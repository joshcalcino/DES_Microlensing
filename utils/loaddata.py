import numpy as np
import os
import matplotlib.pyplot as plt


def load_csv_data(filepath):
    """ We want to load the DES data for each stellar object and save that information in a dictonary.
    Then, we will have a dictionary with the keys as object IDs, with each key referencing another dictionary which has
    keys as the bands we are observing."""

    if os.path.exists(filepath):
        print("Path exits, loading in csv data now..")
        csv_data = np.genfromtxt(fname=filepath, skip_header=1, dtype=None, delimiter=',')

    else:
        print("Error loading csv data. Path does not exist.")
        csv_data = None
        exit()

    tmp_key = csv_data[0][0]
    tmp_dict = {}

    object_dict = {}

    for i in range(0, len(csv_data)):
        tmp_row = csv_data[i]

        print(tmp_row)

        # if tmp_row[3] == 99 or tmp_row[4] == 99:

        if tmp_row[0] == tmp_key:
            print(tmp_row[0])
            print(tmp_row[1])
            if tmp_row[1]+'_mag_obs' in tmp_dict:
                tmp_dict[tmp_row[1] + '_mag_obs'] = np.append(tmp_dict[tmp_row[1] + '_mag_obs'], tmp_row[3])
                tmp_dict[tmp_row[1] + '_mag_err'] = np.append(tmp_dict[tmp_row[1] + '_mag_err'], tmp_row[4])
                tmp_dict[tmp_row[1] + '_obs_time'] = np.append(tmp_dict[tmp_row[1] + '_obs_time'], tmp_row[2])

            else:
                tmp_dict[tmp_row[1] + '_mag_obs'] = np.array([tmp_row[3]])
                tmp_dict[tmp_row[1] + '_mag_err'] = np.array([tmp_row[4]])
                tmp_dict[tmp_row[1] + '_obs_time'] = np.array([tmp_row[2]])

            # print(tmp_dict)

        else:
            object_dict[str(last_tmp_row[0])] = tmp_dict

            tmp_dict = {}
            tmp_key = tmp_row[0]

            tmp_dict[tmp_row[1] + '_mag_obs'] = np.array([tmp_row[3]])
            tmp_dict[tmp_row[1] + '_mag_err'] = np.array([tmp_row[4]])
            tmp_dict[tmp_row[1] + '_obs_time'] = np.array([tmp_row[2]])

        last_tmp_row = tmp_row

    return object_dict

#
# data = load_csv_data('../data/high_nepochs_objects.csv')
#
# key_list = data.keys()
#
# test_dict = data[key_list[300]]
#
# des_filter_arr = ['g', 'r', 'i', 'z', 'Y']
#
# for des_filter in des_filter_arr:
#     plt.errorbar(test_dict[des_filter + '_obs_time'], test_dict[des_filter + '_mag_obs'],
#                  yerr=test_dict[des_filter + '_mag_err'], label=des_filter, fmt='o')
#
# plt.gca().invert_yaxis()  # So that decreasing y-axis is up, consistent with smaller magnitudes being brighter
#
# plt.legend()
# plt.show()

#
# for key in data:
#     print(key)
#
# print("The length of the data dict is: ", len(data))
# print(data[1])