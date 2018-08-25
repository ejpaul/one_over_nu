LIBSTELL_OBJ_FILES = ezcdf.o ezcdf_GenPut.o ezcdf_inqvar.o stel_constants.o v3_utilities.o \
	ezcdf_GenGet.o ezcdf_attrib.o ezcdf_opncls.o handle_err.o read_boozer_mod.o read_wout_mod.o safe_open_mod.o stel_kinds.o vsvd0.o \
	vparams.o mgrid_mod.o vmec_input.o system_mod.o parse_extension.o getcarg.o vmec_getenv.o \
    ezspline.o ezspline_init.o ezspline_obj.o ezspline_type.o ezspline_setup.o \
    akherm1.o akherm2.o akherm3.o ezspline_2netcdf.o splinck.o r8mktricub.o \
    r8mkspline.o r8mkintrp3d.o r8mkintrp2d.o r8mkbicub.o r8genxpkg.o genxpkg.o \
    r8util_bcherm1.o r8util_bcherm2.o r8util_bcherm3.o r8splinck.o r8cspline.o \
    r8akherm1.o r8akherm2.o r8akherm3.o r8v_spline.o mktricub.o mkspline.o \
    mkintrp3d.o mkintrp2d.o cspline.o ibc_ck.o mkbicub.o util_bcherm1.o \
    util_bcherm2.o util_bcherm3.o v_spline.o ezspline_interp.o \
    vecbicub.o vecherm1.o vecherm2.o vecherm3.o evbicub.o evintrp2d.o \
    evintrp3d.o evspline.o xlookup.o vectricub.o vecspline.o vecpc3.o \
    vecpc2.o vecpc1.o vecintrp3d.o vecintrp2d.o vecin3d_argchk.o vecin2d_argchk.o \
    r8vecbicub.o r8vecherm1.o r8vecherm2.o r8vecherm3.o r8vecintrp2d.o r8vecintrp3d.o \
    r8vecpc1.o r8vecpc2.o r8vecpc3.o r8vecspline.o r8vectricub.o r8evintrp2d.o \
    r8evintrp3d.o r8pc1ev.o r8pc2ev.o r8pc3ev.o r8xlookup.o  \
    r8herm3ev.o r8herm2ev.o r8herm1ev.o r8zonfind.o r8gridtricub.o \
    r8gridpc1.o r8gridpc2.o r8gridpc3.o evtricub.o gridbicub.o r8gridherm2.o r8gridherm3.o \
    r8gridintrp2d.o r8gridintrp3d.o gridpc2.o gridpc3.o gridtricub.o herm1ev.o \
    herm2ev.o herm3ev.o pc1ev.o pc2ev.o pc3ev.o r8evbicub.o r8evspline.o \
    r8evtricub.o r8gridbicub.o r8gridherm2.o r8gridherm3.o r8gridintrp2d.o r8gridintrp3d.o \
    gridherm2.o gridherm3.o gridintrp2d.o gridintrp3d.o zonfind.o ezspline_modulo.o FZERO.o 

LIBSTELL_TARGET = mini_libstell.a

.PHONY: all clean

all: $(LIBSTELL_TARGET)

%.o: %.f90
	$(FC) $(EXTRA_COMPILE_FLAGS) -DNETCDF -c $<

%.o: %.f
	$(FC) $(EXTRA_COMPILE_FLAGS) -DNETCDF -c $<

%.o: %.F
	$(FC) $(EXTRA_COMPILE_FLAGS) -DNETCDF -c $<

$(LIBSTELL_TARGET): $(LIBSTELL_OBJ_FILES)
	ar rcs $(LIBSTELL_TARGET) $(LIBSTELL_OBJ_FILES)

clean:
	rm -f *.o *.mod *.MOD *~ $(LIBSTELL_TARGET)

# Dependencies:
FZERO.o:
akherm1.o:
akherm2.o:
akherm3.o:
ezcdf.o: ezcdf_GenPut.o ezcdf_GenGet.o ezcdf_attrib.o ezcdf_opncls.o
ezcdf_GenPut.o: ezcdf_opncls.o ezcdf_inqvar.o
ezcdf_inqvar.o:
stel_constants.o: stel_kinds.o
v3_utilities.o: stel_kinds.o
ezcdf_GenGet.o: ezcdf_opncls.o ezcdf_inqvar.o
ezcdf_attrib.o: ezcdf_inqvar.o
ezcdf_opncls.o:
ezspline.o: ezspline_obj.o
ezspline_init.o: ezspline.o ezspline_setup.o
ezspline_obj.o: ezspline_type.o
ezspline_setup.o:
ezspline_type.o:
ezspline_2netcdf.o:
handle_err.o:
read_boozer_mod.o: stel_kinds.o safe_open_mod.o stel_constants.o ezcdf.o
read_wout_mod.o: stel_kinds.o safe_open_mod.o stel_constants.o v3_utilities.o vsvd0.o ezcdf.o mgrid_mod.o vmec_input.o
safe_open_mod.o:
stel_kinds.o:
vsvd0.o: stel_kinds.o
mgrid_mod.o: stel_kinds.o vmec_input.o vsvd0.o system_mod.o safe_open_mod.o ezcdf.o
vparams.o: stel_kinds.o stel_constants.o
vmec_input.o: vparams.o vsvd0.o
system_mod.o:
parse_extension.o:
getcarg.o:
vmec_getenv.o:
splinck.o:
r8mktricub.o: 
r8mkspline.o:
r8mkintrp3d.o:
r8mkintrp2d.o:
r8mkbicub.o:
r8genxpkg.o:
genxpkg.o:
r8util_bcherm1.o:
r8util_bcherm2.o:
r8util_bcherm3.o:
r8splinck.o:
r8cspline.o:
r8akherm1.o:
r8akherm2.o:
r8akherm3.o:
r8v_spline.o:
mktricub.o:
mkspline.o:
mkintrp2d.o:
mkintrp3d.o:
cspline.o:
ibc_ck.o:
mkbicub.o:
util_bcherm1.o:
util_bcherm2.o:
util_bcherm3.o:
v_spline.o:
ezspline_interp.o:
vecherm1.o:
vecbicub.o:
vecherm2.o:
vecherm3.o:
evbicub.o:
evintrp2d.o:
evintrp3d.o:
evspline.o:
xlookup.o:
vectricub.o:
vecspline.o:
vecpc3.o:
vecpc2.o:
vecpc1.o:
vecintrp3d.o:
vecintrp2d.o:
vecin3d_argchk.o:
vecin2d_argchk.o:
r8vecbicub.o:
r8vecherm1.o:
r8vecherm2.o:
r8vecherm3.o:
r8vecintrp2d.o:
r8vecintrp3d.o:
r8vecpc1.o:
r8vecpc2.o:
r8vecpc3.o:
r8vecspline.o:
r8vectricub.o:
r8pc1ev.o:
r8pc2ev.o:
r8pc3ev.o:
r8xlookup.o:
r8herm3xyz.o:
r8zonfind.o:
r8gridtricub.o:
r8gridpc1.o:
r8gridpc2.o:
r8gridpc3.o:
evtricub.o:
gridbicub.o:
r8gridherm2.o:
r8gridintrp2d.o:
r8gridintrp3d.o:
gridpc2.o:
gridpc3.o:
gridtricub.o:
herm1ev.o:
herm2ev.o:
herm3ev.o:
pc1ev.o:
pc2ev.o:
pc3ev.o:
r8evbicub.o:
r8evspline.o:
r8evtricub.o:
r8gridbicub.o:
r8gridherm2.o:
r8gridherm3.o:
r8gridintrp2d.o:
r8gridintrp3d.o:
gridherm2.o:
gridherm3.o:
gridintrp2d.o:
gridintrp3d.o:
zonfind.o:
ezspline_modulo.o:
