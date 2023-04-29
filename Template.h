#ifndef TEMPLATE_H
#define TEMPLATE_H

#include "constants.h"

template<typename type_pos>
void generateElementpos(type_pos* Elementpos){
	Elementpos[0] = (NUMELEMENTS - 1) / 2 * ELEMENTSPACING;
	for (int i = 1; i < NUMELEMENTS; i++) {
		Elementpos[i] = Elementpos[i-1] + ELEMENTSPACING;

	}
}

template<typename type_din, typename type_dout>
type_dout datatypeConverter(type_din data_i){
	type_dout data_o = (type_dout)data_i;
	return data_o;
}



template<typename type_angle, typename type_dout>
type_dout CalculateSin(type_angle angle)
{
    type_dout output = hls::sinf(angle);
    return output;
}

template<typename type_angle, typename type_dout>
type_dout CalculateCos(type_angle angle)
{
    type_dout output = hls::cosf(angle);
    return output;
}



template<typename din_fix, typename d_int_index, typename din_angle>
class DelayCalculator {
public:
    din_fix computeDelayTime(d_int_index elementIndex, din_angle azimuthAngle, din_angle elevationAngle) {
        din_fix sinAzimuthangle = CalculateSin(azimuthAngle);
        din_fix cosElevationangle = CalculateCos(elevationAngle);
        din_fix delayTime = ((din_fix)((elementIndex) * ELEMENTSPACING / SPEEDOFLIGHT)) * sinAzimuthangle * cosElevationangle;
        return delayTime;
    }

private:
    din_fix CalculateSin(din_angle angle) {
        din_fix functionInput = angle;
        din_fix output = hls::sinf(functionInput);
        return output;
    }

    din_fix CalculateCos(din_angle angle) {
        din_fix functionInput = angle;
        din_fix output = hls::cosf(functionInput);
        return output;
    }


};

















#endif

