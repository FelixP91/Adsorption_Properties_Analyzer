import pandas as pd
from Bio.PDB import PDBParser
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.DSSP import DSSP
from LHP_Gils_2022 import ProteinPatch
import os

pd.set_option('display.max_rows', 10000)
pd.set_option('display.max_columns', 10000)
pd.set_option('display.width', 10000)

# path to the file, containing AA_seqeunces (in fasta-format as downloaded from uniprot) of enzymes. Can be a .txt file.
fasta_file = 'C:\\Users\\panisf91\\Seqeunces_TYRs_wetland.txt'

# path to the folder containing pdb-files of enzymes. Title of pdb files has to be in the format: anything-Uniprot_ID-anything.pdb
# (when pdb files are downloaded from Alphafold, the titles are in the correct format)
pdb_files = 'C:\\Users\\panisf91\\wetland_TYRs_alphafold_pdb_files'

csv_file_name = 'soil_TYRs_csv_12_12_23.txt'  # name of the generated output csv file. The file will be putputted in the currentvworking direvctory.

pHs = [5, 7, 9]  # list with pH values, at which the overall charge and the surface charge will be calculated


class SeqeunceAnalyzer():
    '''This class provides tool to analyze amino cid sequences.
    Input files:
    fasta file containing names and amno acid seqeunces for all analyzed enzymes
    pdb files (from AlphaFold) of all analyzed enzymes.

    Output file: csv file containing standard and optional features of enzymes.
    -standard feaures: UniProt ID, amino acid sequence, Phylum, is-Streptomyces,
    -optional features: LHP (largest hydrophobic patch), LCP (largest charged patch), RHSA (relative hydrophobic surface area), RCSA (relative charged surface area),
    IEP (isoelectric point), charge_at_pHs (total net charge at selected pHs), SSF (secondary structure fraction of loops, helices, sheets), surface charge at pH values.

    optional features can be calculated selectively by calling the correspinding method. By calling complete_analysis() all features are calculated and saved as an output file.
    fasta_file=path to the fast file
    pdb_file_folder= path to folder containing .pdb files. This folder should contain ONLY pdb files.
    output_path= path and name of hte output file. default = 'output_file.txt' in the current working directory
    pH_values= a list of pH values, at which the total charge and the surface charge will be calculated. default = pH 5,7,9.
    '''
    def __init__(self, fasta_file_, pdb_file_folder, output_path='\\'.join([os.getcwd(), 'output_file.txt']), pH_values=[5,7,9]):
        self.fasta_file = fasta_file_
        self.pdb_file_folder = pdb_file_folder
        self.output_path = output_path
        self.pH_values = pH_values
        self.sequences = []
        self.IDs = []

        for enzyme in SeqIO.parse(self.fasta_file, 'fasta'):
            self.sequences.append(enzyme.seq)
            self.IDs.append(enzyme.id)

        self.output_df = pd.DataFrame({'IDs_preliminary': self.IDs, 'sequences': self.sequences})
        self.output_df['split_id'] = self.output_df.IDs_preliminary.apply(lambda x: x.split('_'))
        self.output_df['ID'] = self.output_df.split_id.apply(lambda x: x[0])
        self.output_df['phylum'] = self.output_df.split_id.apply(lambda x: x[1])
        self.output_df['Seq'] = self.output_df.sequences.apply(lambda x: ''.join(x))
        self.output_df['Streptomyces'] = self.output_df.split_id.apply(lambda x: self.is_Streptomyces(x))

        self.output_df.drop(['IDs_preliminary', 'sequences', 'split_id'], axis=1, inplace=True)

        self.LHPs = []
        self.LCPs = []
        self.RHSAs = []
        self.RCSAs = []

    def complete_analysis(self):
        '''Performs a complete analysis of th einput data and creates and output file.'''
        self.charge_at_pHs()
        self.SSF()
        self.IEP()
        self.calculate_patches()
        self.calculate_surface_areas()
        self.calculate_surface_charge()
        self.output_file()


    def is_Streptomyces(self,x):
        '''Returns 1 if an amino acid sequence downloaded from the UniProt Databank belongs to the genus Streptomyces
        and returns 0 if it belongs to any other genus'''
        if len(x) == 3:
            if x[2] == 'Strepto':
                return 1
            else:
                return 0
        else:
            return 0

    def calculate_RHSA(self,x):  # calculate relative hydrophobic surface area
        '''Calculates the surface area ratio, that is covered by hydrophobic amino acids by dividing
        the surface area covered by hydrophobic amino acids by the total surface area of the protein.
        ['A', 'P', 'F', 'I', 'L', 'M', 'V', 'W', 'G'] are defined as hydrophobic amino acids'''
        TSA = x.RSA_A2.sum()
        HSA = x.RSA_A2[x[1].isin(['A', 'P', 'F', 'I', 'L', 'M', 'V', 'W', 'G'])].sum()
        return HSA / TSA

    def calculate_RCSA(self,x):  # calculate relative charged surface area
        '''Calculates the surface area ratio, that is covered by charged amino acids by dividing
                the surface area covered by charged amino acids by the total surface area of the protein.
                ['R', 'D', 'E', 'H', 'K'] are defined as charged amino acids'''
        TSA = x.RSA_A2.sum()
        HSA = x.RSA_A2[x[1].isin(['R', 'D', 'E', 'H', 'K'])].sum()
        return HSA / TSA

    def charge_at_pHs(self):
        '''Calculates the charge of an amino acid sequence at certain pH values and adds them to a DataFrame in the column [f'charge_pH_{pH}']. The default pH values are pH 5,7,9.'''
        for pH in self.pH_values:
            self.output_df[f'charge_pH_{pH}'] = self.output_df.Seq.apply(lambda x: ProteinAnalysis(x).charge_at_pH(pH))

    def SSF(self):
        '''Calculates the secondary structure fraction of an amino acid seqeunce. For each amino acid seqeunce in a DataFrame it
        adds a Column with the ratio of amin acids featured in helices ['SSF_H'], turns ['SSF_T'], and sheets ['SSF_S']'''
        self.output_df['SSF_H'] = self.output_df.Seq.apply(lambda x: ProteinAnalysis(x).secondary_structure_fraction()[0])
        self.output_df['SSF_T'] = self.output_df.Seq.apply(lambda x: ProteinAnalysis(x).secondary_structure_fraction()[1])
        self.output_df['SSF_S'] = self.output_df.Seq.apply(lambda x: ProteinAnalysis(x).secondary_structure_fraction()[2])

    def IEP(self):
        '''Calculates the isoelectric points of amino acid sequences and adds them to a DataFrame in the column ['IEP']'''
        self.output_df['IEP'] = self.output_df.Seq.apply(lambda x: ProteinAnalysis(x).isoelectric_point())

    def calculate_patches(self): #eventuell noch nicht richtig?
        '''Calculates the area (in A^2) of the largest connected hydrophobic and charged patch in th esurface of proteins.
        hydrophobic amino acids: ['A', 'V', 'L', 'I', 'P', 'F', 'M','W']
        charged amino acids: ['R', 'D', 'E', 'H', 'K']. The data is added as new column to a dataframe.
        [LHP] = largest hydrophibic patch.
        [LCP] = largest charged patch.'''
        self.Uniprot_IDs = []
        for root, dirs, files in os.walk(self.pdb_file_folder):
            for file in files:
                file_path_ = os.path.join(root, file)
                # Check if it's a file (not a subdirectory) before processing
                if os.path.isfile(file_path_):
                    Uniprot_ID_ = file_path_.split('\\')[-1].split('-')[1]

                    hydrophobic_ = ProteinPatch(Uniprot_ID_, file_path_, ['A', 'V', 'L', 'I', 'P', 'F', 'M','W'])  # calculate largest hydrophobic patch
                    charged_ = ProteinPatch(Uniprot_ID_, file_path_, ['R', 'D', 'E', 'H', 'K'])  # calculate largest charged patch

                    self.Uniprot_IDs.append(Uniprot_ID_)
                    self.LHPs.append(hydrophobic_.plot_largest_patches())
                    self.LCPs.append(charged_.plot_largest_patches())
        # adding the areas of teh LHP and LCP for each protein to output df. Only the intersection of enzymes provided via a fasta file and as pdf files will be present in the output file.
        self.df_patches = pd.DataFrame({'ID':self.Uniprot_IDs, 'LHP':self.LHPs, 'LCP':self.LCPs})
        self.output_df = pd.merge(self.output_df, self.df_patches, on='ID', how='inner')


    def calculate_surface_areas(self):
        '''Adds columns to a dataframe containing the surface are ratios covered by hydrophobic [RHSA]
        and charged [RCSA] amino acids.'''
        self.Uniprot_IDs_RSA = []
        #accessible surface areas of amino acids
        self.ASA = {'A': 121, 'R': 265, 'N': 187, 'D': 187, 'C': 148, 'E': 214, 'Q': 214, 'G': 97, 'H': 216,
                       'I': 195, 'L': 191, 'K': 230, 'M': 203, 'F': 228, 'P': 154, 'S': 143, 'T': 163, 'W': 264, 'Y': 255, 'V': 165}
        for root, dirs, files in os.walk(self.pdb_file_folder):
            for file in files:
                file_path_ = os.path.join(root, file)
                # Check if it's a file (not a subdirectory) before processing
                if os.path.isfile(file_path_):
                    Uniprot_ID_ = file_path_.split('\\')[-1].split('-')[1]
                    self.p = PDBParser()
                    self.structure = self.p.get_structure(Uniprot_ID_, file_path_)
                    self.model = self.structure[0]
                    #This calculates several parameters for each aminoa cid, including the type of the amino acid(1-letter code) and the accessible surface area ratio.
                    #The data can be read as a dataframe. Column 1 = amino acid (1-letter-code), column 3 = accessible surface area ratio.
                    self.dssp = DSSP(self.model, file_path_)
                    self.df_temp = pd.DataFrame(self.dssp)[[1,3]]
                    #This calculates the accessible surface area by multiplyinf the total surface area with the accessible surface area ratio for each amino acid.
                    self.df_temp['RSA_A2'] = self.df_temp.apply(lambda x: self.ASA[x[1]] * x[3], axis=1)
                    self.RHSAs.append(self.calculate_RHSA(self.df_temp))
                    self.RCSAs.append(self.calculate_RCSA(self.df_temp))
                    self.Uniprot_IDs_RSA.append(Uniprot_ID_)
                    self.df_temp = pd.DataFrame({'ID':self.Uniprot_IDs_RSA,'RHSA':self.RHSAs, 'RCSA':self.RCSAs})
        # adding the RHSA and RCSA values for each protein to output df. Only the intersection of enzymes provided via a fasta file and as pdf files will be present in the output file.
        self.output_df = pd.merge(self.output_df, self.df_temp,on='ID', how='inner')

    def calculate_surface_charge(self):
        '''First, the partial charges of charged amino acids at the defined pH values are calculated. Then, the sum the charges of all surface-exposed amino acids is calculated.
        If more than 20% of the total surface of an amino acid are surface exposed it is classified as a surface-exposed amino acid.'''
        self.pKa_values = {'R': 12.48, 'D': 3.65, 'C': 8.18, 'E': 4.25, 'H': 6, 'K': 10.53, 'Y': 10.07}
        for pH in self.pH_values:
            self.Uniprot_IDs = []
            self.s_charges = []

            #Calculating the partial charges of charged amino acids at specific pH values.
            self.charge_at_pH_dict = {}
            for i in self.pKa_values.keys():
                if i in list('DCEY'):
                    charged_fraction = 10 ** (pH - self.pKa_values[i])
                    self.charge_at_pH_dict[i] = -charged_fraction / (1 + charged_fraction)
                if i in list('RHK'):
                    charged_fraction = 1 / 10 ** (pH - self.pKa_values[i])
                    self.charge_at_pH_dict[i] = charged_fraction / (1 + charged_fraction)

            #loopin over all pdb files and identifying surface exposed amino acids for each enzyme.
            for root, dirs, files in os.walk(self.pdb_file_folder):
                for file in files:
                    self.file_path_ = os.path.join(root, file)
                    # Check if it's a file (not a subdirectory) before processing
                    if os.path.isfile(self.file_path_):
                        self.Uniprot_ID_ = self.file_path_.split('\\')[-1].split('-')[1]
                        self.Uniprot_IDs.append(self.Uniprot_ID_)
                        self.p = PDBParser()
                        self.structure = self.p.get_structure(self.Uniprot_ID_, self.file_path_)
                        self.model = self.structure[0]
                        # This calculates several parameters for each aminoa cid, including the type of the amino acid(1-letter code) and the accessible surface area ratio.
                        # The data can be read as a dataframe. Column 1 = amino acid (1-letter-code), column 3 = accessible surface area ratio.
                        self.dssp = DSSP(self.model, self.file_path_)
                        self.df_temp = pd.DataFrame(self.dssp)[[1, 3]]
                        #Filtering for charged amino acids. If >= 20% of the total surface of an amino acid are surface exposed, the amino acid is classified as a surface exposed amino acid.
                        self.df_temp['surface_charge'] = self.df_temp.apply(lambda x: self.charge_at_pH_dict[x[1]] if x[1] in self.charge_at_pH_dict.keys() and x[3] > 0.2 else 0, axis=1)
                        #Summing the partial charges of all surface exposed charged amino acids to calculate the net toal charge of the protein surface.
                        self.s_charges.append(self.df_temp.surface_charge.sum())
            #adding the surface charges for each protein to output df. Only the intersection of enzymes provided via a fasta file and as pdf files will be present in the output file.
            self.df_temp = pd.DataFrame({'ID': self.Uniprot_IDs, f'S_charges_pH_{pH}': self.s_charges})
            self.output_df = pd.merge(self.output_df, self.df_temp, on='ID', how='inner')

    def output_file(self):
        '''Creates an outpu file. By default, the file "output_file.txt" is created in the cwd. output_df is saved to the file. It contains Amino acid seqeunces of proteins and
        features previously selected.'''
        self.output_df.to_csv(path_or_buf=self.output_path, sep=',')

if __name__ == '__main__':
    x = SeqeunceAnalyzer(fasta_file, pdb_files)
    x.complete_analysis()
    print(x.output_df)

