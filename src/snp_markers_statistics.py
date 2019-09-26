import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.append('/home/ariel/clients/my-tools/src/my_plot_tools')
#import plot_tools

#def show_allele_frequency(df, my_title):
#    sampels_num = df.shape[1] - 2
#    count1 = (df == 1).sum(axis=1)
#    count2 = (df == 2).sum(axis=1)
#    count_all = count1 + count2
#    pr_of_total = count_all / sampels_num
#    maf_freq = pd.concat([count1 / count_all, count2 / count_all], axis=1).min(axis=1)

 #   fig, axs = plt.subplots(1, 2, constrained_layout=True)
 #   fig.suptitle(my_title, fontsize=16)

  #  plot_tools.plot_hist_on_axis(axs[0], pr_of_total, [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], '1 or 2 fraction', 'fraction', 'count')
  #  plot_tools.plot_hist_on_axis(axs[1], maf_freq, None, '1 or 2 fraction', 'fraction', 'count')

   # plt.show()


def compute_informative_markers(df, max_major_allele_freq, min_samples_presence):

    [count1, count2] = get_alleles_count(df)
    count_all = count1 + count2
    informative = (1 * ((count1 / count_all) < max_major_allele_freq)
                   + 2 * ((count2 / count_all) < max_major_allele_freq)) \
                  * (count_all > min_samples_presence)
    return informative


def get_alleles_count(df):

    df1 = df.loc[:, df.columns[2:]]
    count1 = (df1 == 1).sum(axis=1)
    count2 = (df1 == 2).sum(axis=1)
    return [count1, count2]





