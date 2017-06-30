import pandas as pd
# Make the graphs a bit prettier, and bigger
pd.set_option('display.width', 5000) 
pd.set_option('display.max_columns', 60) 
#%%
import getpass
user = getpass.getuser()
path = '/home/{:}/Documents/thesis/data/LiteratureReview/TableWLE/'.format(user)
dbFile = 'wleLiterature.csv'
db = pd.read_csv(path+dbFile,sep='\t')
print(db[:0])
#%%

is_2d = db['2D'] == 1
is_sym = db['Aerofoil'].str.startswith('NACA00',na=False)
is_postStall = db['Post-Stall'] == '1'

selection = db[is_2d & ~is_sym & is_postStall]


print(selection[['Author','Year']])
