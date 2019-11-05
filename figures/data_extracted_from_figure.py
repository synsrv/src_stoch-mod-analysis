

x00 = (85.5, 320.5)

y08 = (85.5, 79.)

x12 = (259.5, 320.5)


d_x0 = (85.5, 18.5)
d_x0_up = (85.5, 18.5)
d_x0_dw = (85.5, 18.5)


d_x4 = (143.5, 200.5)
d_x4_up = (143.5, 197.5)
d_x4_dw = (143.5, 205.5)
    

d_x8 = (201.5, 253.5)
d_x8_up = (201.5, 249.5)
d_x8_dw = (201.5, 257.5)

d_x12 = (259.5, 275.5)
d_x12_up = (259.5, 270.5)
d_x12_dw = (259.5, 279.5)

d_x16 = (317.5, 283.5)
d_x16_up = (317.5, 276.5)
d_x16_dw = (317.5, 289.5)


dps_pix = [d_x0, d_x4, d_x8, d_x12, d_x16]
dps_up = [d_x0_up, d_x4_up, d_x8_up, d_x12_up, d_x16_up]
dps_dw = [d_x0_dw, d_x4_dw, d_x8_dw, d_x12_dw, d_x16_dw]

xvals_data = [(dp[0]-x00[0])/(x12[0]-x00[0])*12 for dp in dps_pix]
yvals_data = [(dp[1]-x00[1])/(y08[1]-x00[1])*0.8 for dp in dps_pix]

yvals_data=[round(y,3) for y in yvals_data]


yerrs = [((dps_up[i][1]-dps_dw[i][1])/2.)/(y08[1]-x00[1])*0.8 for i in range(5)]


if __name__ == "__main__":

    print(xvals_data, yvals_data)
