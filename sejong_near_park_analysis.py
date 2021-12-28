import os 
import geopandas as gpd
from tqdm import tqdm 
import haversine as hs
from haversine import Unit
import pandas as pd 

''' Data Load ''' 
data_path = os.path.join(os.getcwd(), 'data')
gdf_31 = gpd.read_file(os.path.join(data_path, '31.세종시_법정경계(읍면동).geojson'))
gdf_23 = gpd.read_file(os.path.join(data_path, '23.세종시_도로명주소_건물.geojson'))
park_gdf = gpd.read_file(os.path.join(data_path, 'TL_SPOT_PARK.shp'), encoding='euc-kr')

'''
<전처리 1>
생활권에 해당하는 면폴리곤만 남기기
- gdf_31은 세종시 전체를 포함하고 있는 데이터이므로, 1~6생활권 폴리곤만 추출해야함 
'''
happy_city_gdf = gdf_31.loc[(gdf_31['EMD_KOR_NM'] == '집현동') | (gdf_31['EMD_KOR_NM'] == '반곡동') | (gdf_31['EMD_KOR_NM'] == '소담동') | (gdf_31['EMD_KOR_NM'] == '보람동') | 
                            (gdf_31['EMD_KOR_NM'] == '대평동') | (gdf_31['EMD_KOR_NM'] == '가람동') | (gdf_31['EMD_KOR_NM'] == '한솔동') | (gdf_31['EMD_KOR_NM'] == '새롬동') | (gdf_31['EMD_KOR_NM'] == '나성동') |
                            (gdf_31['EMD_KOR_NM'] == '다정동') | (gdf_31['EMD_KOR_NM'] == '어진동') | (gdf_31['EMD_KOR_NM'] == '종촌동') | (gdf_31['EMD_KOR_NM'] == '고운동') | (gdf_31['EMD_KOR_NM'] == '도담동') |
                            (gdf_31['EMD_KOR_NM'] == '아름동') | (gdf_31['EMD_KOR_NM'] == '해밀동') | (gdf_31['EMD_KOR_NM'] == '산울동') | (gdf_31['EMD_KOR_NM'] == '합강동')]

# r(region) 
r1 = ['고운동', '아름동', '종촌동', '도담동', '어진동'] # 1생활권
r2 = ['다정동', '새롬동', '한솔동', '나성동', '가람동'] # 2생활권
r3 = ['대평동', '보람동', '소담동'] # 3생활권
r4 = ['반곡동', '집현동'] # 4생활권
r5 = ['합강동'] # 5생활권
r6 = ['산울동', '해밀동'] # 6생활권 

rs = []
for i in range(len(happy_city_gdf)): 
    cur_emd_cd = happy_city_gdf.iloc[i]['EMD_CD']
    cur_emd_nm = happy_city_gdf.iloc[i]['EMD_KOR_NM']
    
    if cur_emd_nm in r1: 
        rs.append('1생활권')
    if cur_emd_nm in r2: 
        rs.append('2생활권')
    if cur_emd_nm in r3: 
        rs.append('3생활권')
    if cur_emd_nm in r4: 
        rs.append('4생활권')
    if cur_emd_nm in r5: 
        rs.append('5생활권')
    if cur_emd_nm in r6: 
        rs.append('6생활권')
happy_city_gdf['생활권'] = rs
rs.append('중앙녹지공원')




'''
<전처리 2>
생활권 별 건물 나누기
'''
gdf_23 = gdf_23.sample(frac=0.8).reset_index(drop=True) # 건물 개수가 너무 많아서 샘플링 했다. 

# 1생활권
gdf_23_1 = gdf_23.loc[(gdf_23['EMD_CD'] == '111') | (gdf_23['EMD_CD'] == '112') | (gdf_23['EMD_CD'] == '113') | (gdf_23['EMD_CD'] == '114') | (gdf_23['EMD_CD'] == '110')]

# 2생활권
gdf_23_2 = gdf_23.loc[(gdf_23['EMD_CD'] == '106') | (gdf_23['EMD_CD'] == '107') | (gdf_23['EMD_CD'] == '108') | (gdf_23['EMD_CD'] == '109') | (gdf_23['EMD_CD'] == '105')]

# 3생활권
gdf_23_3 = gdf_23.loc[(gdf_23['EMD_CD'] == '102') | (gdf_23['EMD_CD'] == '103') | (gdf_23['EMD_CD'] == '104')]

# 4생활권
gdf_23_4 = gdf_23.loc[(gdf_23['EMD_CD'] == '101') | (gdf_23['EMD_CD'] == '118')]

# 5생활권
gdf_23_5 = gdf_23.loc[(gdf_23['EMD_CD'] == '117')]

# 6생활권 
gdf_23_6 = gdf_23.loc[(gdf_23['EMD_CD'] == '115') | (gdf_23['EMD_CD'] == '116')]






# 좌표계 통일 (5179 --> 4326), 세종시에 해당하는 녹지공간만 추출 
park_gdf.set_crs(epsg=5179, inplace=True, allow_override=True)
park_gdff = park_gdf.loc[(park_gdf['SIG_CD'] == '36110')] # 세종시의 시군구코드는 36110 
park_gdff_ = park_gdff.to_crs(epsg=4326)


''' 
각 생활권 별로 ‘집에서’ 녹지공간까지의 접근성이 동일한가? -- happy_city_park, happy_city_park_dong
'''
print('각 생활권 별로 ‘집에서’ 녹지공간까지의 접근성이 동일한가? -- happy_city_park \n')

def assign_dong_value(dictionary, df_dong, attr): 
    
    df_dong[attr] = ''
    for i in range(len(df_dong)):
        cur_dong_nm = df_dong.iloc[i]['EMD_NM']
        
        for k, v in dictionary.items(): 
            if k == cur_dong_nm : 
                df_dong.at[i, attr] = v
                
    return df_dong

def calculate_dist_to_park(gdf, park_gdff_): 
    
    bldg_2_park = []
    
    for i in tqdm(range(len(gdf))):
        # 해당 gdf(동)내 모든 빌딩에 대하여 중심점을 구함
        cur_bldg_coord_x = gdf.iloc[i]['geometry'].centroid.xy[0][0]
        cur_bldg_coord_y = gdf.iloc[i]['geometry'].centroid.xy[1][0] 
        
        park_dist_ls  = []
        min_d = 9999
        for j in range(len(park_gdff_)): 
            # 세종시 내 모든 녹지공원의 중심점을 구함 
            cur_park_x = park_gdff_.iloc[j]['geometry'].centroid.xy[0][0]
            cur_park_y = park_gdff_.iloc[j]['geometry'].centroid.xy[1][0] 
            
            meter_dist = hs.haversine((cur_park_y, cur_park_x), 
                                      (cur_bldg_coord_y, cur_bldg_coord_x), 
                                      unit=Unit.METERS)
            
            if (meter_dist/1000.0) < 5: 
                polygon_type = park_gdff_.iloc[j]['geometry'].type
                
                # park_gdf의 폴리곤이 두 종류로 이루어져 있어, 종류에 따라 다른 계산법을 적용해야함 
                if polygon_type == 'MultiPolygon': 
                    for q in range(len(list(park_gdff_.iloc[j]['geometry']))):
                        polys = list(park_gdff_.iloc[j]['geometry'])
                        polyy = polys[q]
                        polyy_crd_length = len(polyy.exterior.coords)     
                        
                        for mc in range(polyy_crd_length):
                            ex_x = polyy.exterior.coords.xy[0][mc]
                            ex_y = polyy.exterior.coords.xy[1][mc]
                            
                            dist1 = hs.haversine((ex_y, ex_x), (cur_bldg_coord_y, cur_bldg_coord_x), unit=Unit.METERS)
                            if dist1 < min_d: 
                                min_d = dist1
                        park_dist_ls.append(min_d)
                            
                if polygon_type == 'Polygon':
                    simple_polygon = park_gdff_.iloc[j]['geometry'].simplify(0.05)
                    for c in range(len(simple_polygon.exterior.coords.xy[0])):
                
                        ex_xx = simple_polygon.exterior.coords.xy[0][c]
                        ex_yy = simple_polygon.exterior.coords.xy[1][c]
                
                        dist2 = hs.haversine((ex_yy, ex_xx), (cur_bldg_coord_y, cur_bldg_coord_x), unit=Unit.METERS)
                    
                        if dist2 < min_d: 
                            min_d = dist2
                    park_dist_ls.append(min_d)
                        
            else:       
                park_dist_ls.append(min_d)
                
        bldg_2_park.append(min(park_dist_ls))
    
    return bldg_2_park

def find_dong_nm(gdf, happy_city_gdf):
    
    gdf['DONG_NM'] = ''
    
    for i in range(len(gdf)):
        cur_apt_poly = gdf.iloc[i]['geometry']
        
        if cur_apt_poly == None : 
            continue 
        
        cur_apt_centroid = cur_apt_poly[0].centroid 

        
        for j in range(len(happy_city_gdf)):
            cur_dong_poly = happy_city_gdf.iloc[j]['geometry'][0]
            cur_dong_nm = happy_city_gdf.iloc[j]['EMD_KOR_NM']
            
            if cur_dong_poly.contains(cur_apt_centroid): 
                gdf.at[i, 'DONG_NM'] = cur_dong_nm

    
    return gdf 

            

bldg_2_park_gdf_23_1 = calculate_dist_to_park(gdf_23_1, park_gdff_)
gdf_23_1 = gdf_23_1.reset_index()
gdf_23_1 = gdf_23_1.drop(['index'], axis=1)
gdf_23_1 = find_dong_nm(gdf_23_1, happy_city_gdf)
gdf_23_1['최근린녹지거리'] = bldg_2_park_gdf_23_1

bldg_2_park_gdf_23_2 = calculate_dist_to_park(gdf_23_2, park_gdff_)
gdf_23_2 = gdf_23_2.reset_index()
gdf_23_2 = gdf_23_2.drop(['index'], axis=1)
gdf_23_2 = find_dong_nm(gdf_23_2, happy_city_gdf)
gdf_23_2['최근린녹지거리'] = bldg_2_park_gdf_23_2

bldg_2_park_gdf_23_3 = calculate_dist_to_park(gdf_23_3, park_gdff_)
gdf_23_3 = gdf_23_3.reset_index()
gdf_23_3 = gdf_23_3.drop(['index'], axis=1)
gdf_23_3 = find_dong_nm(gdf_23_3, happy_city_gdf)
gdf_23_3['최근린녹지거리'] = bldg_2_park_gdf_23_3

bldg_2_park_gdf_23_4 = calculate_dist_to_park(gdf_23_4, park_gdff_)
gdf_23_4 = gdf_23_4.reset_index()
gdf_23_4 = gdf_23_4.drop(['index'], axis=1)
gdf_23_4 = find_dong_nm(gdf_23_4, happy_city_gdf)
for i in range(len(gdf_23_4)):
    cur_ = gdf_23_4.iloc[i]['EMD_CD']
    if cur_ == '101': 
        gdf_23_4.at[i, 'DONG_NM'] = '반곡동'
gdf_23_4['최근린녹지거리'] = bldg_2_park_gdf_23_4

bldg_2_park_gdf_23_5 = calculate_dist_to_park(gdf_23_5, park_gdff_)
gdf_23_5 = gdf_23_5.reset_index()
gdf_23_5 = gdf_23_5.drop(['index'], axis=1)
gdf_23_5 = find_dong_nm(gdf_23_5, happy_city_gdf)
gdf_23_5['최근린녹지거리'] = bldg_2_park_gdf_23_5

bldg_2_park_gdf_23_6 = calculate_dist_to_park(gdf_23_6, park_gdff_)
gdf_23_6 = gdf_23_6.reset_index()
gdf_23_6 = gdf_23_6.drop(['index'], axis=1)
gdf_23_6 = find_dong_nm(gdf_23_6, happy_city_gdf)
gdf_23_6['최근린녹지거리'] = bldg_2_park_gdf_23_6


happy_city_park = {'1생활권': gdf_23_1['최근린녹지거리'].sum()/len(gdf_23_1), 
                   '2생활권': gdf_23_2['최근린녹지거리'].sum()/len(gdf_23_2),
                   '3생활권': gdf_23_3['최근린녹지거리'].sum()/len(gdf_23_3),
                   '4생활권': gdf_23_4['최근린녹지거리'].sum()/len(gdf_23_4),
                   '5생활권': gdf_23_5['최근린녹지거리'].sum()/len(gdf_23_5),
                   '6생활권': gdf_23_6['최근린녹지거리'].sum()/len(gdf_23_6)}
    
happy_city_park_dong = {'고운동': gdf_23_1.loc[gdf_23_1['DONG_NM']=='고운동']['최근린녹지거리'].sum()/len(gdf_23_1[gdf_23_1['DONG_NM']=='고운동']), 
                        '아름동': gdf_23_1.loc[gdf_23_1['DONG_NM']=='아름동']['최근린녹지거리'].sum()/len(gdf_23_1[gdf_23_1['DONG_NM']=='아름동']),
                        '종촌동': gdf_23_1.loc[gdf_23_1['DONG_NM']=='종촌동']['최근린녹지거리'].sum()/len(gdf_23_1[gdf_23_1['DONG_NM']=='종촌동']),
                        '도담동': gdf_23_1.loc[gdf_23_1['DONG_NM']=='도담동']['최근린녹지거리'].sum()/len(gdf_23_1[gdf_23_1['DONG_NM']=='도담동']),
                        '어진동': gdf_23_1.loc[gdf_23_1['DONG_NM']=='고운동']['최근린녹지거리'].sum()/len(gdf_23_1[gdf_23_1['DONG_NM']=='어진동']),
                        '다정동': gdf_23_2.loc[gdf_23_2['DONG_NM']=='다정동']['최근린녹지거리'].sum()/len(gdf_23_2[gdf_23_2['DONG_NM']=='다정동']),
                        '새롬동': gdf_23_2.loc[gdf_23_2['DONG_NM']=='새롬동']['최근린녹지거리'].sum()/len(gdf_23_2[gdf_23_2['DONG_NM']=='새롬동']),
                        '한솔동': gdf_23_2.loc[gdf_23_2['DONG_NM']=='한솔동']['최근린녹지거리'].sum()/len(gdf_23_2[gdf_23_2['DONG_NM']=='한솔동']),
                        '나성동': gdf_23_2.loc[gdf_23_2['DONG_NM']=='나성동']['최근린녹지거리'].sum()/len(gdf_23_2[gdf_23_2['DONG_NM']=='나성동']),
                        '가람동': gdf_23_2.loc[gdf_23_2['DONG_NM']=='가람동']['최근린녹지거리'].sum()/len(gdf_23_2[gdf_23_2['DONG_NM']=='가람동']),
                        '대평동': gdf_23_3.loc[gdf_23_3['DONG_NM']=='대평동']['최근린녹지거리'].sum()/len(gdf_23_3[gdf_23_3['DONG_NM']=='대평동']),
                        '보람동': gdf_23_3.loc[gdf_23_3['DONG_NM']=='보람동']['최근린녹지거리'].sum()/len(gdf_23_3[gdf_23_3['DONG_NM']=='보람동']),
                        '소담동': gdf_23_3.loc[gdf_23_3['DONG_NM']=='소담동']['최근린녹지거리'].sum()/len(gdf_23_3[gdf_23_3['DONG_NM']=='소담동']),
                        '반곡동': gdf_23_4.loc[gdf_23_4['DONG_NM']=='반곡동']['최근린녹지거리'].sum()/len(gdf_23_4[gdf_23_4['DONG_NM']=='반곡동']),
                        '집현동': gdf_23_4.loc[gdf_23_4['DONG_NM']=='집현동']['최근린녹지거리'].sum()/len(gdf_23_4[gdf_23_4['DONG_NM']=='집현동']),
                        '합강동': gdf_23_5.loc[gdf_23_5['DONG_NM']=='합강동']['최근린녹지거리'].sum()/len(gdf_23_5[gdf_23_5['DONG_NM']=='합강동']),
                        '산울동': gdf_23_6.loc[gdf_23_6['DONG_NM']=='산울동']['최근린녹지거리'].sum()/len(gdf_23_6[gdf_23_6['DONG_NM']=='산울동']),
                        '해밀동': gdf_23_6.loc[gdf_23_6['DONG_NM']=='해밀동']['최근린녹지거리'].sum()/len(gdf_23_6[gdf_23_6['DONG_NM']=='해밀동'])}    
                  
dong_dict = {'고운동': [0, 0], '아름동': [0, 0], '종촌동': [0, 0], '도담동': [0, 0], '어진동': [0, 0], '다정동': [0, 0],
             '새롬동': [0, 0], '한솔동': [0, 0], '나성동': [0, 0], '가람동': [0, 0], '대평동': [0, 0], '보람동': [0, 0],
             '소담동': [0, 0], '반곡동': [0, 0], '집현동': [0, 0], '합강동': [0, 0], '산울동': [0, 0], '해밀동': [0, 0]}

df_dong = pd.DataFrame(dong_dict.keys(), columns=['EMD_NM'])                                                     
df_dong = assign_dong_value(happy_city_park_dong, df_dong, '최근린녹지까지의평균거리') 
