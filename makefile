.phony: clean veryclean save plot_erg

clean:
	rm -f dump_[0-9]* pmd_[0-9]* *~

veryclean: clean
	rm -f out.* {erg,frc,strs}.pmd

resdir := result_$(shell date "+%y%m%d_%H%M")
save:
	mkdir -p $(resdir)
	cp in.* out.* pmdini pmdfin dump_* $(resdir)/

plot_erg:
	echo "plot 'out.erg' us 1:3 w l t 'total', '' us 1:4 w l t 'kinetic', '' us 1:5 w l t 'potential'" > plot_erg.gp
	gnuplot -e "load 'plot_erg.gp'; pause -1"
