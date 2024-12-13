
import pandas as pd
import numpy as np
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

data_table = pd.read_csv("variants.tsv", sep="\t");

# Soru 1 Frequency sütunundaki eksik değerleri, o varyantın bulunduğu popülasyonın ortalama frekansı ile doldurun.
means_of_populations_frequency = data_table.groupby('Population')['Frequency'].mean()
print(means_of_populations_frequency)

# Boş değerler,tablodaki satırın, popülasyon değerinin ortalamasıyla doldurulur.
data_table['Frequency'] = data_table.apply( lambda row: means_of_populations_frequency[row['Population']] if np.isnan(row['Frequency']) else row['Frequency'], axis=1 )
# Eğer sütundaki NaN değerler de ortalamaya dahil edilecekse:
# NaN değerler 0 ile doldurulur
# data_table['Frequency'] = data_table['Frequency'].fillna(0.0000)
# Ortalama alınır
# means_of_populations_with_zero = data_table.groupby('Population')['Frequency'].mean()
# Boş değerler,tablodaki satırın, popülasyon değerinin ortalamasıyla doldurulur.
# data_table['Frequency'] = data_table.apply(
#     lambda row: means_of_populations_with_zero[row['Population']] if row['Frequency'] == 0.0000 else row['Frequency'],
#     axis=1
# )

# Soru 2 Clinical_Significance sütunundaki eksik değerleri, o popülasyondaki en sık görülen Clinical_Significance değeriyle doldurun.
# Popülasyon'a göre gruplanmış Clinical_Significance değerlerinin en çok tekrar eden değeri bulunur
most_significance_for_population = data_table.groupby('Population')['Clinical_Significance'].apply(
    lambda x:
    x.mode().iloc[0]
    if not x.mode().empty
    else None
)
# tabloda Clinical_Significance değeri boş olan satırlar, hesaplanan, populasyona göre en çok tekrar eden Clinical Significance değeriyle doldurulur.
data_table['Clinical_Significance'] = data_table.apply(
    lambda row:
    most_significance_for_population[row['Population']]
    if pd.isna(row['Clinical_Significance'])
    else row['Clinical_Significance'],
    axis=1
)
# Soru 3 Her popülasyonun ortalama frekans değerini hesaplayın. Ve bu ortalama değerin üzerinde ve eşit frekansa sahip olan varyantları filtreleyin.
# Yalnızca Pathogenic olarak işaretlenmiş varyantları içeren bir alt küme oluşturun.

# Yukarıda hesaplamış olduğumuz, popülasyonun frekans ortalamalarına göre, tablodaki satırların ortalama değerlerine bakılır ve yüksek ve eşit ise filtrelenir.
filtered_mean = data_table[data_table.apply(
    lambda row: row["Frequency"] >= means_of_populations_frequency[row["Population"]]
    if pd.notna(row["Frequency"])
    else False,
    axis=1)
]
#sadece pathogenic varyanların alt kümesi oluşturulur.
pathogenic_variants = filtered_mean[filtered_mean["Clinical_Significance"] == "Pathogenic"]
sub_set_pathogenic = set(pathogenic_variants['Variant'])
# Soru 4 Filtrelenmiş varyantlar üzerinde, bir Venn diyagramı ile varyantların popülasyonlar arasında ortaklık bilgisinin görselleştirilmesini sağlayın.
# tablodaki bütün popülasyon sütunundaki tüm değerler bulunur.
####################################v
# Dinamik olarak hesaplama
#all_populations = data_table["Population"].unique()
#tablodaki populasyon değerine göre varyantların alt kümesi oluşturulur
#variant_sets = {pop: set(data_table[data_table["Population"] == pop]["Variant"]) for pop in all_populations}
# alt kümeler venn şemasında gösterilmesi için listeye atılır
#sets_list = list(variant_sets.values())
#if len(sets_list) == 3:
#    venn = venn3(
#        sets_list,
#        set_labels=all_populations
#    )
#    plt.title("Venn Diagram of Variants by Population")
#    plt.show()
#else:
#    print("Venn Diagram not draw. Population count: {len(all_populations)}")
####################################v
# Bilinen Populasyon değerlerine göre hesaplama
#
european_variants = set(data_table[data_table["Population"] == "European"]["Variant"])
african_variants = set(data_table[data_table["Population"] == "African"]["Variant"])
asian_variants = set(data_table[data_table["Population"] == "Asian"]["Variant"])
venn = venn3(
    [european_variants, african_variants, asian_variants],
    ("European", "African", "Asian")
)
plt.title("Venn Diagram with Variants")
plt.show()
