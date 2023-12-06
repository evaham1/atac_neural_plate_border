# set up work dir
work_dir = '/home/output/NF-downstream_analysis/Processing/ss8/SEACELLS_ATAC_WF/2_SEACells_computation/exported_data'
os.listdir(work_dir)

# not sure what this is for yet
# if not os.path.exists(os.path.join(work_dir, 'model')):
#     os.makedirs(os.path.join(work_dir, 'model'))

# read in the summarised peak counts across SEACells
mat = pd.read_csv(os.path.join(work_dir, "Summarised_by_metacells_counts.csv"))

# reformat matrix
mat = mat.set_index('Unnamed: 0') # turn first column into rownames
mat = mat.iloc[1: , :] # delete empty row
mat = mat.transpose() # transpose table so seacells are columns and peaks are rows

# create cisTopic oboject from the matrix
cisTopic_obj = create_cistopic_object(mat)