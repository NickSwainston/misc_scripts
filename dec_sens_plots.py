import find_pulsar_in_obs as fpio
import numpy as np
import matplotlib.pyplot as plt

def calc_powers(centrefreq=150.):
    sweet_dec_range = [-82.8,-71.4,-63.1,-55.,-47.5,-40.4,-33.5,-26.7,-19.9,-13.,-5.9,1.6,9.7,18.6,29.4,44.8]
    sweet_delays_range= [[0,0,0,0,7,7,7,7,14,14,14,14,21,21,21,21],\
                         [0,0,0,0,6,6,6,6,12,12,12,12,18,18,18,18],\
                         [0,0,0,0,5,5,5,5,10,10,10,10,15,15,15,15],\
                         [0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12],\
                         [0,0,0,0,3,3,3,3,6,6,6,6,9,9,9,9],\
                         [0,0,0,0,2,2,2,2,4,4,4,4,6,6,6,6],\
                         [0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3],\
                         [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],\
                         [3,3,3,3,2,2,2,2,1,1,1,1,0,0,0,0],\
                         [6,6,6,6,4,4,4,4,2,2,2,2,0,0,0,0],\
                         [9,9,9,9,6,6,6,6,3,3,3,3,0,0,0,0],\
                         [12,12,12,12,8,8,8,8,4,4,4,4,0,0,0,0],\
                         [15,15,15,15,10,10,10,10,5,5,5,5,0,0,0,0],\
                         [18,18,18,18,12,12,12,12,6,6,6,6,0,0,0,0],\
                         [21,21,21,21,14,14,14,14,7,7,7,7,0,0,0,0],\
                         [24,24,24,24,16,16,16,16,8,8,8,8,0,0,0,0]]

    #Going to work out how many pointings are needed
    #setting up some metadata requirements
    time = 4800 #one hour 20 min
    channels = range(int(centrefreq/1.28)-12, int(centrefreq/1.28)+12)
    minfreq = float(min(channels))
    maxfreq = float(max(channels))

    obsid = '1117624530'
    ra = 180.
    dec_calc = []
    ra_calc = []
    for i in range(-89,89,1):
            dec_calc.append(i)
            ra_calc.append(ra)
    
    powers_list = []
    for delays in sweet_delays_range:
        ra_sex, deg_sex = fpio.deg2sex(ra, dec_calc[0])
        cord = [obsid, str(ra_sex), str(deg_sex), 1, delays, centrefreq, channels]
        #powout=get_beam_power(cord, zip(RA_FWHM_calc,Dec_FWHM_calc), dt=600)
        names_ra_dec = np.column_stack((['source']*len(ra_calc), ra_calc, dec_calc))
        powout = fpio.get_beam_power_over_time(cord, names_ra_dec, dt=600, degrees=True)
        temp_power = []
        for p in powout:
            temp_power.append(p[0][0])
        powers_list.append(temp_power)
        #print(temp_power)
    
    return powers_list

if __name__ == "__main__":
    for freq in [120, 150, 200]:
        powers_list = calc_powers(freq)
        sweet_dec_range = [-82.8,-71.4,-63.1,-55.,-47.5,-40.4,-33.5,-26.7,-19.9,-13.,-5.9,1.6,9.7,18.6,29.4,44.8]
        plt.figure(figsize=(15,10))
        for pi, p in enumerate(powers_list):
            if sweet_dec_range[pi] < -30.:
                plt.plot(range(-89,89,1), p, linestyle=":", label=sweet_dec_range[pi])
            else:
                plt.plot(range(-89,89,1), p, label=sweet_dec_range[pi])
        plt.legend()
        plt.xlabel("Declination (degrees)")
        plt.ylabel("Zenith normalised power")
        plt.title("Meridian powers for observations at {} MHz".format(freq))
        plt.xlim(-90, 65)
        plt.savefig("dec_sens_{}.png".format(freq))

    for freq in [120, 150, 200]:
        powers_list = calc_powers(freq)
        sweet_dec_range = [-82.8,-71.4,-63.1,-55.,-47.5,-40.4,-33.5,-26.7,-19.9,-13.,-5.9,1.6,9.7,18.6,29.4,44.8]
        plt.figure(figsize=(15,10))
        for pi, p in enumerate(powers_list):
            #check if they're worth using above 30 deg
            if max(p[120:]) < 0.05:
                continue
            
            if sweet_dec_range[pi] < -30.:
                plt.plot(range(-89,89,1), p, linestyle=":", label=sweet_dec_range[pi])
            else:
                plt.plot(range(-89,89,1), p, label=sweet_dec_range[pi])
        plt.legend()
        plt.xlabel("Declination (degrees)")
        plt.ylabel("Zenith normalised power")
        plt.title("Meridian powers for observations at {} MHz".format(freq))
        plt.xlim(30, 65)
        plt.ylim(0.,0.3)
        plt.savefig("dec_sens_zoom_{}.png".format(freq))


   
