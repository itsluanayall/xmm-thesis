import os
from config import CONFIG
from astropy.table import Table, vstack

def combine_obs_tables(products_dir, output_name):

    #Initialize final table and its path
    final_table = Table()
    final_table_path = os.path.join(products_dir, output_name)

    #Go through the observations
    for obsid in os.listdir(products_dir):
        if obsid.startswith('0'):

            print(f'Processing obs {obsid}.')

            #Read table of single observation
            table_path = os.path.join(products_dir, f"{obsid}", f"{obsid}_table.fits")
            print(table_path)
            table = Table.read(table_path, format='fits')

            #Append to final table
            final_table = vstack([final_table, table])
            final_table.write(final_table_path, format='fits', overwrite=True)
            print(f'Appended observation {obsid}!')

        else:
            continue

    return final_table
        
if __name__ == "__main__":

    #target_dir = CONFIG['target_dir']
    target_dir = "/home/luana/Desktop/Magistrale/Thesis/Markarian421"
    products_dir = os.path.join(target_dir, "Products", "RGS_Spectra")
    spectra_table = combine_obs_tables(products_dir=products_dir, output_name="spectra_table.fits")
