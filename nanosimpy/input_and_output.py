import tifffile as tifffile
import numpy as np
import uuid

def save_to_tiff(data_store, fwhm, path_and_file_name):
    """Saves out data as a tiff file, do a given directory.
    --inputs--
    data_store:         is the dictionary containing the intensity traces.
    fwhm:               is the array index of the fwhm being exported.
    path_and_filename:   the destination and filename of the tiff-file being saved (folder/exp.tif)
    --outputs--
    Saves tiff to directory.

    """
    #Creates matrix big enough to hold entire carpet.
    
    matrix_im = np.zeros((data_store.__len__(), data_store[0]['trace'][fwhm].shape[0]-4))
    for row_idx in data_store:
        #Populates matrix with values.
        matrix_im[row_idx, :] = data_store[row_idx]['trace'][fwhm][:matrix_im.shape[1]]
    with tifffile.TiffWriter(path_and_file_name, bigtiff=True) as tif:
        tif.save(matrix_im.T.astype(np.float64))
    return
def save_as_csv(time_series, data, path_and_file_name):
    """Save out individual data track
    --inputs--
    time_series:        The time points associated with data.
    data:               The correlation data.
    path_and_file_name:  the destination path and filename.
    --outputs--
    Saves the file as a .csv file conmpatible with FoCuS-point and FoCuS-scan fitting software.
    """

    parent_name = "test"
    parent_uqid = uuid.uuid4()

    file_obj = open(path_and_file_name, 'w')
    file_obj.write('version,'+str(2)+'\n')
    file_obj.write('numOfCH,'+str(1)+'\n')
    file_obj.write('type, scan\n')
    file_obj.write('ch_type,'+str(0)+'\n')

    file_obj.write('carpet pos,'+str(0)+'\n')
    file_obj.write('parent_name,'+str(parent_name)+'\n')
    file_obj.write('parent_uqid,'+str(parent_uqid)+'\n')
    file_obj.write('parent_filename,'+str(path_and_file_name)+'\n')

    file_obj.write('pc, 0\n')
    file_obj.write('Time (ns), CH0 Auto-Correlation\n')
    for time_step in range(0, time_series.shape[0]):
        file_obj.write(str(float(time_series[time_step]))+','+str(data[time_step])+ '\n')
    file_obj.write('end\n')

    file_obj.close()