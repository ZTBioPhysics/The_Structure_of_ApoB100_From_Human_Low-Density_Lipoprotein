#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created: April 2024

@author: Zachary Berndsen, ChatGPT
Description: 
    This script can be run to create all the PhoX crosslinking plots and tables that appear in the paper
"""

import xlink_analysis_functions as xlf
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option('display.max_columns', None)  # None means show all columns

##############################################################################################################################################################################

apob_seq = 'MDPPRPALLALLALPALLLLLLAGARAEEEMLENVSLVCPKDATRFKHLRKYTYNYEAESSSGVPGTADSRSATRINCKVELEVPQLCSFILKTSQCTLKEVYGFNPEGKALLKKTKNSEEFAAAMSRYELKLAIPEGKQVFLYPEKDEPTYILNIKRGIISALLVPPETEEAKQVLFLDTVYGNCSTHFTVKTRKGNVATEISTERDLGQCDRFKPIRTGISPLALIKGMTRPLSTLISSSQSCQYTLDAKRKHVAEAICKEQHLFLPFSYKNKYGMVAQVTQTLKLEDTPKINSRFFGEGTKKMGLAFESTKSTSPPKQAEAVLKTLQELKKLTISEQNIQRANLFNKLVTELRGLSDEAVTSLLPQLIEVSSPITLQALVQCGQPQCSTHILQWLKRVHANPLLIDVVTYLVALIPEPSAQQLREIFNMARDQRSRATLYALSHAVNNYHKTNPTGTQELLDIANYLMEQIQDDCTGDEDYTYLILRVIGNMGQTMEQLTPELKSSILKCVQSTKPSLMIQKAAIQALRKMEPKDKDQEVLLQTFLDDASPGDKRLAAYLMLMRSPSQADINKIVQILPWEQNEQVKNFVASHIANILNSEELDIQDLKKLVKEALKESQLPTVMDFRKFSRNYQLYKSVSLPSLDPASAKIEGNLIFDPNNYLPKESMLKTTLTAFGFASADLIEIGLEGKGFEPTLEALFGKQGFFPDSVNKALYWVNGQVPDGVSKVLVDHFGYTKDDKHEQDMVNGIMLSVEKLIKDLKSKEVPEARAYLRILGEELGFASLHDLQLLGKLLLMGARTLQGIPQMIGEVIRKGSKNDFFLHYIFMENAFELPTGAGLQLQISSSGVIAPGAKAGVKLEVANMQAELVAKPSVSVEFVTNMGIIIPDFARSGVQMNTNFFHESGLEAHVALKAGKLKFIIPSPKRPVKLLSGGNTLHLVSTTKTEVIPPLIENRQSWSVCKQVFPGLNYCTSGAYSNASSTDSASYYPLTGDTRLELELRPTGEIEQYSVSATYELQREDRALVDTLKFVTQAEGAKQTEATMTFKYNRQSMTLSSEVQIPDFDVDLGTILRVNDESTEGKTSYRLTLDIQNKKITEVALMGHLSCDTKEERKIKGVISIPRLQAEARSEILAHWSPAKLLLQMDSSATAYGSTVSKRVAWHYDEEKIEFEWNTGTNVDTKKMTSNFPVDLSDYPKSLHMYANRLLDHRVPQTDMTFRHVGSKLIVAMSSWLQKASGSLPYTQTLQDHLNSLKEFNLQNMGLPDFHIPENLFLKSDGRVKYTLNKNSLKIEIPLPFGGKSSRDLKMLETVRTPALHFKSVGFHLPSREFQVPTFTIPKLYQLQVPLLGVLDLSTNVYSNLYNWSASYSGGNTSTDHFSLRARYHMKADSVVDLLSYNVQGSGETTYDHKNTFTLSCDGSLRHKFLDSNIKFSHVEKLGNNPVSKGLLIFDASSSWGPQMSASVHLDSKKKQHLFVKEVKIDGQFRVSSFYAKGTYGLSCQRDPNTGRLNGESNLRFNSSYLQGTNQITGRYEDGTLSLTSTSDLQSGIIKNTASLKYENYELTLKSDTNGKYKNFATSNKMDMTFSKQNALLRSEYQADYESLRFFSLLSGSLNSHGLELNADILGTDKINSGAHKATLRIGQDGISTSATTNLKCSLLVLENELNAELGLSGASMKLTTNGRFREHNAKFSLDGKAALTELSLGSAYQAMILGVDSKNIFNFKVSQEGLKLSNDMMGSYAEMKFDHTNSLNIAGLSLDFSSKLDNIYSSDKFYKQTVNLQLQPYSLVTTLNSDLKYNALDLTNNGKLRLEPLKLHVAGNLKGAYQNNEIKHIYAISSAALSASYKADTVAKVQGVEFSHRLNTDIAGLASAIDMSTNYNSDSLHFSNVFRSVMAPFTMTIDAHTNGNGKLALWGEHTGQLYSKFLLKAEPLAFTFSHDYKGSTSHHLVSRKSISAALEHKVSALLTPAEQTGTWKLKTQFNNNEYSQDLDAYNTKDKIGVELTGRTLADLTLLDSPIKVPLLLSEPINIIDALEMRDAVEKPQEFTIVAFVKYDKNQDVHSINLPFFETLQEYFERNRQTIIVVLENVQRNLKHINIDQFVRKYRAALGKLPQQANDYLNSFNWERQVSHAKEKLTALTKKYRITENDIQIALDDAKINFNEKLSQLQTYMIQFDQYIKDSYDLHDLKIAIANIIDEIIEKLKSLDEHYHIRVNLVKTIHDLHLFIENIDFNKSGSSTASWIQNVDTKYQIRIQIQEKLQQLKRHIQNIDIQHLAGKLKQHIEAIDVRVLLDQLGTTISFERINDILEHVKHFVINLIGDFEVAEKINAFRAKVHELIERYEVDQQIQVLMDKLVELAHQYKLKETIQKLSNVLQQVKIKDYFEKLVGFIDDAVKKLNELSFKTFIEDVNKFLDMLIKKLKSFDYHQFVDETNDKIREVTQRLNGEIQALELPQKAEALKLFLEETKATVAVYLESLQDTKITLIINWLQEALSSASLAHMKAKFRETLEDTRDRMYQMDIQQELQRYLSLVGQVYSTLVTYISDWWTLAAKNLTDFAEQYSIQDWAKRMKALVEQGFTVPEIKTILGTMPAFEVSLQALQKATFQTPDFIVPLTDLRIPSVQINFKDLKNIKIPSRFSTPEFTILNTFHIPSFTIDFVEMKVKIIRTIDQMLNSELQWPVPDIYLRDLKVEDIPLARITLPDFRLPEIAIPEFIIPTLNLNDFQVPDLHIPEFQLPHISHTIEVPTFGKLYSILKIQSPLFTLDANADIGNGTTSANEAGIAASITAKGESKLEVLNFDFQANAQLSNPKINPLALKESVKFSSKYLRTEHGSEMLFFGNAIEGKSNTVASLHTEKNTLELSNGVIVKINNQLTLDSNTKYFHKLNIPKLDFSSQADLRNEIKTLLKAGHIAWTSSGKGSWKWACPRFSDEGTHESQISFTIEGPLTSFGLSNKINSKHLRVNQNLVYESGSLNFSKLEIQSQVDSQHVGHSVLTAKGMALFGEGKAEFTGRHDAHLNGKVIGTLKNSLFFSAQPFEITASTNNEGNLKVRFPLRLTGKIDFLNNYALFLSPSAQQASWQVSARFNQYKYNQNFSAGNNENIMEAHVGINGEANLDFLNIPLTIPEMRLPYTIITTPPLKDFSLWEKTGLKEFLKTTKQSFDLSVKAQYKKNKHRHSITNPLAVLCEFISQSIKSFDRHFEKNRNNALDFVTKSYNETKIKFDKYKAEKSHDELPRTFQIPGYTVPVVNVEVSPFTIEMSAFGYVFPKAVSMPSFSILGSDVRVPSYTLILPSLELPVLHVPRNLKLSLPDFKELCTISHIFIPAMGNITYDFSFKSSVITLNTNAELFNQSDIVAHLLSSSSSVIDALQYKLEGTTRLTRKRGLKLATALSLSNKFVEGSHNSTVSLTTKNMEVSVATTTKAQIPILRMNFKQELNGNTKSKPTVSSSMEFKYDFNSSMLYSTAKGAVDHKLSLESLTSYFSIESSTKGDVKGSVLSREYSGTIASEANTYLNSKSTRSSVKLQGTSKIDDIWNLEVKENFAGEATLQRIYSLWEHSTKNHLQLEGLFFTNGEHTSKATLELSPWQMSALVQVHASQPSSFHDFPDLGQEVALNANTKNQKIRWKNEVRIHSGSFQSQVELSNDQEKAHLDIAGSLEGHLRFLKNIILPVYDKSLWDFLKLDVTTSIGRRQHLRVSTAFVYTKNPNGYSFSIPVKVLADKFIIPGLKLNDLNSVLVMPTFHVPFTDLQVPSCKLDFREIQIYKKLRTSSFALNLPTLPEVKFPEVDVLTKYSQPEDSLIPFFEITVPESQLTVSQFTLPKSVSDGIAALDLNAVANKIADFELPTIIVPEQTIEIPSIKFSVPAGIVIPSFQALTARFEVDSPVYNATWSASLKNKADYVETVLDSTCSSTVQFLEYELNVLGTHKIEDGTLASKTKGTFAHRDFSAEYEEDGKYEGLQEWEGKAHLNIKSPAFTDLHLRYQKDKKGISTSAASPAVGTVGMDMDEDDDFSKWNFYYSPQSSPDKKLTIFKTELRVRESDEETQIKVNWEEEAASGLLTSLKDNVPKATGVLYDYVNKYHWEHTGLTLREVSSKLRRNLQNNAEWVYQGAIRQIDDIDVRFQKAASGTTGTYQEWKDKAQNLYQELLTQEGQASFQGLKDNVFDGLVRVTQEFHMKVKHLIDSLIDFLNFPRFQFPGKPGIYTREELCTMFIREVGTVLSQVYSKVHNGSEILFSYFQDLVITLPFELRKHKLIDVISMYRELLKDLSKEAQEVFKAIQSLKTTEVLRNLQDLLQFIFQLIEDNIKQLKEMKFTYLINYIQDEINTIFSDYIPYVFKLLKENLCLNLHKFNEFIQNELQEASQELQQIHQYIMALREEYFDPSIVGWTVKYYELEEKIVSLIKNLLVALKDFHSEYIVSASNFTSQLSSQVEQFLHRNIQEYLSILTDPDGKGKEKIAELSATAQEIIKSQAIATKKIISDYHQQFRYKLQDFSDQLSDYYEKFIAESKRLIDLSIQNYHTFLIYITELLKKLQSTTVMNPYMKLAPGELTIIL'

##############################################################################################################################################################################

# start from annotated csv file

xlinks_all = 'all_common_xlinks_small_and_large.csv'
# xlinks_small = 'all_xlinks_small.csv'
# xlinks_large = 'all_xlinks_large.csv'

df = xlf.read_csv_with_header(xlinks_all)

##############################################################################################################################################################################
# make summaryt plots

# make new data frames that are filtered by spectral count
spectral_count_threshold = 10
filtered_df = xlf.filter_by_spectral_count(df, spectral_count_threshold)
spectral_count_threshold = 20
filtered_df2 = xlf.filter_by_spectral_count(df, spectral_count_threshold)

# CA Dist vs spectral count
plt.figure(figsize=(10, 6))
plt.scatter(df['CA Distance'],df['Spectral Count'], alpha=0.7)
plt.yscale('log')
plt.grid(True)
plt.axvline(x=20, color='r', linestyle='--', linewidth=2)
plt.axhline(y=20, color='r', linestyle='--', linewidth=2)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.xlabel('C-alpha Distance (Å)', fontsize=16)
plt.ylabel('Spectral Count', fontsize=16)
plt.title('C-alpha Distance vs Spectral Count', fontsize=16)
plt.show()

# Seq Dist vs CA Dist
plt.figure(figsize=(10, 6))
plt.scatter(df['Sequence Distance'], df['CA Distance'], alpha=0.7, label="All Data")
plt.scatter(filtered_df2['Sequence Distance'], filtered_df2['CA Distance'], alpha=0.7, color='r', label="Spectral Count >= 20")
plt.grid(True)
plt.axhline(y=26, color='red', linestyle='--', linewidth=2)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.xlabel('Sequence Distance', fontsize=16)
plt.ylabel('C-alpha Distance (Å)', fontsize=16)
plt.title('Sequence Distance vs C-alpha Distance', fontsize=16)
plt.legend(fontsize=14)
plt.show()

# Seq dist vs spec count
plt.figure(figsize=(10, 6))
plt.scatter(df['Sequence Distance'], df['Spectral Count'], alpha=0.7)
plt.yscale('log')
plt.axhline(y=20, color='r', linestyle='--', linewidth=2)
plt.grid(True)
plt.tick_params(axis='both', which='major', labelsize=14)
plt.xlabel('Sequence Distance', fontsize=16)
plt.ylabel('Spectral Count', fontsize=16)
plt.title('Sequence Distance vs Spectral Count', fontsize=16)
plt.legend(fontsize=14)
plt.show()

# cumulative percentage (changed function to min CA Distance)
plt.figure(figsize=(10, 6))  # Create a figure before calling the function
xlf.plot_ca_distance_cumulative_percentage(df, label="All Data")
xlf.plot_ca_distance_cumulative_percentage(filtered_df, label="Spectral Count >= 10")
xlf.plot_ca_distance_cumulative_percentage(filtered_df2, label="Spectral Count >= 20")
plt.axvline(x=26, color='r', linestyle='--')
plt.tick_params(axis='both', which='major', labelsize=14)
plt.xlabel('C-alpha Distance (Å)', fontsize=16)
plt.ylabel('Cumulative Percentage (%)', fontsize=16)
plt.legend(fontsize=14)
plt.title('Cumulative Percentage of Crosslinks <= C-alpha Distance', fontsize=16)
plt.show() 

##############################################################################################################################################################################

# Summarize by domain
summary_df = xlf.domain_crosslink_stats(filtered_df2, apob_seq)
print(summary_df)

# save new df as csv file
output_file = 'all_common_summary_by_domain_sc20.csv'
xlf.save_dataframe_to_csv(summary_df, output_file)

# Summarize by domain association
summary_df = xlf.summarize_by_domain_association(filtered_df2)
print(summary_df)

# save new df as csv file
output_file = 'all_common_summary_by_domaina_ssociation_sc20.csv'
xlf.save_dataframe_to_csv(summary_df, output_file)

##############################################################################################################################################################################

# filter by domain (can be single domain or domain-domain interaction)
selection_string = 'insert 9 to insert 9'  # or simply omit this argument in the call
domains_to_exclude = ['insert 6'] # exlude these domains

filtered_df3 = xlf.filter_by_domain(
    df=filtered_df2,
    domain_string=selection_string,
    exclude_intra_domain=False,
    exclude=False,
    domains_to_exclude=domains_to_exclude
)

print(filtered_df3)

##############################################################################################################################################################################

# print out ChimeraX distance marker commands
mod_number = 1
formatted_list = xlf.format_chimera_dist_selection(filtered_df3, mod_number)
print(formatted_list)

# Copy and past this command in chimeraX to set the display parameters for the crosslinks
# distance style color green radius 0.5 dashes 0

##############################################################################################################################################################################
# other filters

# filter by residue
residue = 4207
filtered_df3 = xlf.filter_by_residue(filtered_df2, residue)
print(filtered_df3)

# filter by Ca-distance
comparison_type = 'less' # can specify less or greater
ca_distance_threshold = 100
filtered_df3 = xlf.filter_by_ca_distance(filtered_df2, ca_distance_threshold, comparison_type)
print(filtered_df3)

##############################################################################################################################################################################
# calculate column statistics

average_ca_distance, stdev_ca_distance = xlf.calculate_average_and_stdev_of_column(filtered_df3, 'CA Distance')
if average_ca_distance is not None:
    print(f"Average CA Distance: {average_ca_distance}, Standard Deviation: {stdev_ca_distance}")

##############################################################################################################################################################################