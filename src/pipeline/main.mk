.PHONY: MAKE_DIRS all


SCRIPTS = $(wildcard ../scripts/*.sh)
TARGET_SCRIPTS = $(patsubst ../scripts/%, ${TARGET_DIR}/%, $(SCRIPTS))
PLGD_PM = $(wildcard ./pipeline/Plgd/*.pm)
TARGET_PLGD_PM = $(patsubst ./pipeline/Plgd/%, ${TARGET_DIR}/Plgd/%, $(PLGD_PM))

all: ${TARGET_DIR}/mecat.pl ${TARGET_DIR}/mecat.sh \
               ${TARGET_PLGD_PM} \
	       ${TARGET_SCRIPTS}

${TARGET_DIR}/mecat.pl: pipeline/mecat.pl
	@if [ ! -e ${TARGET_DIR} ] ; then mkdir -p ${TARGET_DIR}/ ; fi
	cp -pf pipeline/mecat.pl ${TARGET_DIR}/mecat.pl
	chmod +x ${TARGET_DIR}/mecat.pl

${TARGET_DIR}/mecat.sh: pipeline/mecat.sh
	@if [ ! -e ${TARGET_DIR} ] ; then mkdir -p ${TARGET_DIR} ; fi
	cp -pf pipeline/mecat.sh ${TARGET_DIR}/mecat.sh
	chmod +x ${TARGET_DIR}/mecat.sh


$(TARGET_PLGD_PM):${TARGET_DIR}/Plgd/% : ./pipeline/Plgd/% 
	@if [ ! -e ${TARGET_DIR}/Plgd ] ; then mkdir -p ${TARGET_DIR}/Plgd ; fi
	cp -pf  $< $@

$(TARGET_SCRIPTS):${TARGET_DIR}/% : ../scripts/%
	@if [ ! -e ${TARGET_DIR} ] ; then mkdir -p ${TARGET_DIR} ; fi
	cp -pf  $^ $@
	chmod +x $@
	

clean:
	rm -f ${TARGET_DIR}/mecat.pl
	rm -f ${TARGET_DIR}/mecat.sh
	rm -f ${TARGET_SCRIPTS}
	rm -f ${TARGET_PLGD_PM}
	rm -rf ${TARGET_DIR}/Plgd
