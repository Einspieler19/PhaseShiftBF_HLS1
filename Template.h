#ifndef TEMPLATE_H
#define TEMPLATE_H

#include <cstdlib>
#include "constants.h"


template<typename type_data>
struct my_complex_Value{
type_data re;
type_data im;
};

template<typename type_data, int N>
struct my_complex_Array{
type_data re[N];
type_data im[N];
};



template<typename type_pos>
void generateElementpos(type_pos* Elementpos){
	Elementpos[0] = (NUMELEMENTS - 1) / 2 * ELEMENTSPACING;
	for (int i = 1; i < NUMELEMENTS; i++) {
		Elementpos[i] = Elementpos[i-1] + ELEMENTSPACING;

	}
}

template<typename type_din, typename type_dout>
type_dout datatypeConverter(type_din data_i){
	type_dout data_o= 0.0;
	if ((data_i > 1.0e-9)|(data_i < -1.0e-9)){
		data_o = (type_dout)data_i;
	}
	//else{
	//	std::cout<<"no"<< data_i << std::endl;
	//}
	return data_o;
}

template<typename type_din, typename type_dout>
my_complex_Value<type_dout> complexdataConverter(type_din data_i_re, type_din data_i_im){
    my_complex_Value<type_dout> data_o;
    data_o.re = datatypeConverter<type_din, type_dout>(data_i_re);
    data_o.im = datatypeConverter<type_din, type_dout>(data_i_im);
    return data_o;
}


template<typename type_din, typename type_dout, int N>
my_complex_Array<type_dout, N> complexarrayConverter(type_din data_i_re[N], type_din data_i_im[N]){
    my_complex_Array<type_dout, N> data_o;
for (int j = 0; j < N; j++) {
    data_o.re[j] = datatypeConverter<type_din, type_dout>(data_i_re[j]);
    data_o.im[j] = datatypeConverter<type_din, type_dout>(data_i_im[j]);
}
    return data_o;
}


template<typename type_data>
void computeManifoldvector(
		type_data Angle,
		type_data* manifoldRe,
		type_data* manifoldIm)
{
	// 初始化天线阵向量
    type_data Elementpos[NUMELEMENTS];
    generateElementpos<type_data>(Elementpos);

	// 求权重向量
    for (int i = 0; i < NUMELEMENTS; i++) {
    	manifoldRe[i] = hls::cos(2 * M_PI * hls::sinf(Angle) * Elementpos[i] * SPEEDOFLIGHT / CENTERFREQ);
        manifoldIm[i] = hls::sin(2 * M_PI * hls::sinf(Angle) * Elementpos[i] * SPEEDOFLIGHT / CENTERFREQ);
    }
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

