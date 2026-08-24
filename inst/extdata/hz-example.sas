* Excerpted (de-identified) from the public example                          ;
* hazard/examples/hz.death.AVC.sas -- no patient data, only the PROC HAZARD  ;
* block, kept for hzr_translate_sas() documentation.                        ;
%HAZARD(
PROC HAZARD DATA=AVCS CONSERVE OUTHAZ=EXAMPLES.HZDEATH
     QUASI CONDITION=14;
     EVENT DEAD;
     TIME INT_DEAD;
     PARMS MUE=0.2361727 THALF=0.1512095 NU=1.438652 M=1 FIXM
           MUC=0.0005436977;
);
