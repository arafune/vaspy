FC = ifort

TARGET: WaveTrans.x WaveTransSpinor.x WaveTransPlot.x

WaveTrans.x: WaveTrans.f90
	$(FC) -assume byterecl -o WaveTrans.x  WaveTrans.f90

WaveTransSpinor.x: WaveTransSpinor.f90
	$(FC) -assume byterecl  -o WaveTransSpinor.x $<

WaveTransPlot.x: WaveTransPlot.f90
	$(FC) -assume byterecl  -o WaveTransTransplot.x $<
