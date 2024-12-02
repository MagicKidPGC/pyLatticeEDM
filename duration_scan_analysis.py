import io
import zipfile
import numpy as np
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def duration_scan(FILENAME):
    dur = []
    on_counts = []
    off_counts = []
    countsToVms=0.01; # the clock period is 10 us and the ToF yAxis are in the unit of volts
    with zipfile.ZipFile(FILENAME+".zip") as zf:
        for item in zf.namelist():
            if 'average' in item:
                with io.TextIOWrapper(zf.open(item), encoding="utf-8") as f:
                    tree = ET.parse(f)
                    #calling the root element
                    root = tree.getroot()
                    points = root[1]
                    for scanpoint in points:
                        # save the value of the scanparameter
                        dur.append(float(scanpoint[1].text))
                        
                        on_data = scanpoint[2]
                        off_data = scanpoint[3]

                        
                        counts_in_on_shots = []
                        for shot in on_data:
                            ToF = shot[0][0]
                            start = float(ToF[2].text)
                            clockperiod = float(ToF[3].text)

                            tof, counts, bg = [], [], []
                            for i, datapoint in enumerate(ToF[1]):
                                tof_point = start + i*clockperiod
                                
                                if apply_tof_gate:
                                    if (tof_point>=tof_gate_start) & (tof_point<=tof_gate_end):
                                        tof.append(tof_point)
                                        counts.append(float(datapoint.text))
                                    elif (tof_point>=bg_tof_gate_start) & (tof_point<=bg_tof_gate_end):
                                        bg.append(float(datapoint.text))
                                else:
                                    tof.append(tof_point)
                                    counts.append(float(datapoint.text))
                                
                            # plt.plot(tof, counts)
                            # plt.show()

                            counts = np.asarray(counts)
                            c_sum = np.sum(counts)

                            if apply_tof_gate:
                                bg = np.asarray(bg)
                                bg_sum = np.sum(bg)
                                bg_eff=bg_sum/(bg_tof_gate_end-bg_tof_gate_start)*(tof_gate_end-tof_gate_start)
                                counts_in_on_shots.append(c_sum-bg_eff)

                            else:
                                counts_in_on_shots.append(c_sum)
                        counts_in_on_shots = np.asarray(counts_in_on_shots)
                        on_counts.append(np.mean(counts_in_on_shots)*countsToVms)

                        counts_in_off_shots = []
                        for shot in off_data:
                            ToF = shot[0][0]
                            start = float(ToF[2].text)
                            clockperiod = float(ToF[3].text)

                            tof, counts, bg = [], [], []
                            for i, datapoint in enumerate(ToF[1]):
                                tof_point = start + i*clockperiod
                                
                                if apply_tof_gate:
                                    if (tof_point>=tof_gate_start) & (tof_point<=tof_gate_end):
                                        tof.append(tof_point)
                                        counts.append(float(datapoint.text))
                                    elif (tof_point>=bg_tof_gate_start) & (tof_point<=bg_tof_gate_end):
                                        bg.append(float(datapoint.text))
                                else:
                                    tof.append(tof_point)
                                    counts.append(float(datapoint.text))
                                
                            # plt.plot(tof, counts)
                            # plt.show()

                            counts = np.asarray(counts)
                            c_sum = np.sum(counts)

                            if apply_tof_gate:
                                bg = np.asarray(bg)
                                bg_sum = np.sum(bg)
                                bg_eff=bg_sum/(bg_tof_gate_end-bg_tof_gate_start)*(tof_gate_end-tof_gate_start)
                                counts_in_off_shots.append(c_sum-bg_eff)

                            else:
                                counts_in_off_shots.append(c_sum)
                        counts_in_off_shots = np.asarray(counts_in_off_shots)
                        off_counts.append(np.mean(counts_in_off_shots)*countsToVms)


                # print(len(dur), len(on_counts), len(off_counts))

                durarr = np.asarray(dur)
                oncountsarr = np.asarray(on_counts)
                offcountsarr = np.asarray(off_counts)

                return durarr, oncountsarr, offcountsarr


#folder = 'C:\\Users\\u0140952\\Desktop\\CCM\\LatticeEDM\\analysis\\September2024\\2\\'

apply_tof_gate = True
tof_gate_start = 1900
tof_gate_end = 2400
bg_tof_gate_start = 5000
bg_tof_gate_end = 6000


scan1 = '002_DurationScan_V0_001'
duration1, oncounts1, offcounts1 = duration_scan(scan1)
ratio1 = oncounts1 / offcounts1
difference1 = oncounts1 - offcounts1 - np.min(oncounts1 - offcounts1)

print(oncounts1)
print(offcounts1)
print(ratio1)
plt.figure(0)
plt.plot(duration1, oncounts1)
plt.plot(duration1, offcounts1)
plt.xlabel(r'Slowing duration ($\mu$s)')
plt.ylabel('gated ToF area (V.ms)')
plt.legend(["onShots","offShots"])
plt.show()
# plt.plot(duration1, ratio1)
# plt.xlabel(r'Slowing duration ($\mu$s)')
# plt.ylabel('Ratio on_counts / off_counts')
# plt.show()

v0 = 0.932
v1 = 0.065
v2 = 0.0028
v3 = 0.000195
global survival
survival = v0 + v1 + v2 + v3

def func(x, R, bg):
    x = x*1e-6
    return bg + survival**(R*x)

initial_guess = [8e5, 1.2]
params, covariance = curve_fit(func, duration1, ratio1, p0=initial_guess)
R, bg = params
print('Scattering rate: %.2e' % R)
print('Background: %.2e' % bg)
print('Total number of scattered photons:', R*duration1[-1]*1e-6)

# # Generate fitted y data
plt.figure(1)
ratio1_fit = func(duration1, *params)
plt.scatter(duration1, ratio1)
plt.plot(duration1, ratio1_fit)
plt.xlabel(r'Slowing duration ($\mu$s)')
plt.ylabel('Ratio on_counts / off_counts')
plt.xlim(0,100)
plt.ylim(0,2)
plt.show()
print(ratio1)
