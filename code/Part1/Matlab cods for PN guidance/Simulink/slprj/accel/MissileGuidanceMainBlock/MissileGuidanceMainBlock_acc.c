#include "__cf_MissileGuidanceMainBlock.h"
#include <math.h>
#include "MissileGuidanceMainBlock_acc.h"
#include "MissileGuidanceMainBlock_acc_private.h"
#include <stdio.h>
#include "simstruc.h"
#include "fixedpoint.h"
#define CodeFormat S-Function
#define AccDefine1 Accelerator_S-Function
static void mdlOutputs ( SimStruct * S , int_T tid ) { real_T gyxeey3glt ;
real_T * lastU ; real_T c51uyxa2vr ; real_T dixqeke22x ; real_T cpxnpbppov ;
real_T ekr52whu3s ; n3qi1whofz * _rtB ; loikxjbxjg * _rtP ; ew10rzwqr2 *
_rtDW ; _rtDW = ( ( ew10rzwqr2 * ) ssGetRootDWork ( S ) ) ; _rtP = ( (
loikxjbxjg * ) ssGetDefaultParam ( S ) ) ; _rtB = ( ( n3qi1whofz * )
_ssGetBlockIO ( S ) ) ; if ( ssIsSampleHit ( S , 1 , 0 ) ) { c51uyxa2vr =
_rtP -> P_0 ; _rtB -> plqtsmcsxc = _rtP -> P_1 ; } _rtB -> a4sot4u5ht = ( (
f1xhd02yjc * ) ssGetContStates ( S ) ) -> bq3isy3pyn ; _rtB -> pl44ouvwxi = (
( f1xhd02yjc * ) ssGetContStates ( S ) ) -> hnn54tmcbv ; gyxeey3glt = _rtB ->
a4sot4u5ht - _rtB -> pl44ouvwxi ; dixqeke22x = gyxeey3glt * gyxeey3glt ; _rtB
-> oll3pexgfv = ( ( f1xhd02yjc * ) ssGetContStates ( S ) ) -> gopmpownp5 ;
_rtB -> fm25ytexbt = ( ( f1xhd02yjc * ) ssGetContStates ( S ) ) -> nru1t3asy3
; cpxnpbppov = _rtB -> oll3pexgfv - _rtB -> fm25ytexbt ; ekr52whu3s =
cpxnpbppov * cpxnpbppov + dixqeke22x ; if ( ssIsMajorTimeStep ( S ) ) { if (
_rtDW -> mur5fnftlp != 0 ) { ssSetSolverNeedsReset ( S ) ; _rtDW ->
mur5fnftlp = 0 ; } _rtB -> pmhai3cdnu = muDoubleScalarSqrt ( ekr52whu3s ) ; }
else { if ( ekr52whu3s < 0.0 ) { _rtB -> pmhai3cdnu = - muDoubleScalarSqrt (
muDoubleScalarAbs ( ekr52whu3s ) ) ; } else { _rtB -> pmhai3cdnu =
muDoubleScalarSqrt ( ekr52whu3s ) ; } _rtDW -> mur5fnftlp = 1 ; } if ( (
_rtDW -> nln3kpaw1s >= ssGetT ( S ) ) && ( _rtDW -> ooly1cqdaf >= ssGetT ( S
) ) ) { ekr52whu3s = 0.0 ; } else { dixqeke22x = _rtDW -> nln3kpaw1s ; lastU
= & _rtDW -> p2tj0q31lr ; if ( _rtDW -> nln3kpaw1s < _rtDW -> ooly1cqdaf ) {
if ( _rtDW -> ooly1cqdaf < ssGetT ( S ) ) { dixqeke22x = _rtDW -> ooly1cqdaf
; lastU = & _rtDW -> gsdvdewgvg ; } } else { if ( _rtDW -> nln3kpaw1s >=
ssGetT ( S ) ) { dixqeke22x = _rtDW -> ooly1cqdaf ; lastU = & _rtDW ->
gsdvdewgvg ; } } ekr52whu3s = ( _rtB -> pmhai3cdnu - * lastU ) / ( ssGetT ( S
) - dixqeke22x ) ; } ekr52whu3s *= _rtP -> P_6 ; _rtB -> ja2qzpgi0c = (
ekr52whu3s < kesz3x3i1i ( S ) -> fyjmaqcpxz ) ; if ( ssIsSampleHit ( S , 1 ,
0 ) && _rtB -> ja2qzpgi0c ) { ssSetStopRequested ( S , 1 ) ; }
ssCallAccelRunBlock ( S , 0 , 17 , SS_CALL_MDL_OUTPUTS ) ;
ssCallAccelRunBlock ( S , 0 , 18 , SS_CALL_MDL_OUTPUTS ) ;
ssCallAccelRunBlock ( S , 0 , 19 , SS_CALL_MDL_OUTPUTS ) ;
ssCallAccelRunBlock ( S , 0 , 20 , SS_CALL_MDL_OUTPUTS ) ; _rtB -> ohgp1nre5z
= muDoubleScalarAtan ( gyxeey3glt / cpxnpbppov ) ; if ( ( _rtDW -> it4xhky51a
>= ssGetT ( S ) ) && ( _rtDW -> dymjgrmn00 >= ssGetT ( S ) ) ) { cpxnpbppov =
0.0 ; } else { dixqeke22x = _rtDW -> it4xhky51a ; lastU = & _rtDW ->
fcuedkogxl ; if ( _rtDW -> it4xhky51a < _rtDW -> dymjgrmn00 ) { if ( _rtDW ->
dymjgrmn00 < ssGetT ( S ) ) { dixqeke22x = _rtDW -> dymjgrmn00 ; lastU = &
_rtDW -> dpq5ryflf5 ; } } else { if ( _rtDW -> it4xhky51a >= ssGetT ( S ) ) {
dixqeke22x = _rtDW -> dymjgrmn00 ; lastU = & _rtDW -> dpq5ryflf5 ; } }
cpxnpbppov = ( _rtB -> ohgp1nre5z - * lastU ) / ( ssGetT ( S ) - dixqeke22x )
; } if ( ssIsSampleHit ( S , 1 , 0 ) ) { _rtB -> cv54m4abgw = c51uyxa2vr /
_rtB -> plqtsmcsxc ; } gyxeey3glt = ( ( f1xhd02yjc * ) ssGetContStates ( S )
) -> kfubregyho ; if ( ssIsSampleHit ( S , 1 , 0 ) ) { _rtB -> eukpgiz1ov =
_rtP -> P_9 ; _rtB -> msxbp2zh53 = _rtP -> P_10 * _rtB -> plqtsmcsxc ; _rtB
-> llv3ok5rtu = _rtP -> P_11 ; } ekr52whu3s = ekr52whu3s * cpxnpbppov * _rtB
-> llv3ok5rtu ; dixqeke22x = _rtP -> P_12 * ekr52whu3s ; c51uyxa2vr =
muDoubleScalarSin ( _rtB -> ohgp1nre5z + gyxeey3glt ) * _rtB -> plqtsmcsxc /
_rtB -> eukpgiz1ov ; if ( c51uyxa2vr > 1.0 ) { c51uyxa2vr = 1.0 ; } else { if
( c51uyxa2vr < - 1.0 ) { c51uyxa2vr = - 1.0 ; } } cpxnpbppov = _rtB ->
ohgp1nre5z + muDoubleScalarAsin ( c51uyxa2vr ) ; _rtB -> l5zouhgmne =
muDoubleScalarCos ( cpxnpbppov ) * _rtB -> eukpgiz1ov ; if ( _rtDW ->
nn3yzop52g . IcNeedsLoading ) { ( ( f1xhd02yjc * ) ssGetContStates ( S ) ) ->
nfvvmkgt4f = _rtB -> l5zouhgmne ; } _rtB -> afl0c4ktu2 = ( ( f1xhd02yjc * )
ssGetContStates ( S ) ) -> nfvvmkgt4f ; _rtB -> hotk4q5jzw =
muDoubleScalarSin ( cpxnpbppov ) * _rtB -> eukpgiz1ov ; if ( _rtDW ->
k4peu225td . IcNeedsLoading ) { ( ( f1xhd02yjc * ) ssGetContStates ( S ) ) ->
iyebu4ash1 = _rtB -> hotk4q5jzw ; } _rtB -> jcktelybde = ( ( f1xhd02yjc * )
ssGetContStates ( S ) ) -> iyebu4ash1 ; _rtB -> hlobdoowm0 =
muDoubleScalarSin ( gyxeey3glt ) * _rtB -> plqtsmcsxc ; _rtB -> mivzl3g1kg =
muDoubleScalarCos ( gyxeey3glt ) * _rtB -> msxbp2zh53 ; _rtB -> ptld5lrfqw =
muDoubleScalarSin ( _rtB -> ohgp1nre5z ) * dixqeke22x ; _rtB -> nydpnr35ew =
ekr52whu3s * muDoubleScalarCos ( _rtB -> ohgp1nre5z ) ; UNUSED_PARAMETER (
tid ) ; }
#define MDL_UPDATE
static void mdlUpdate ( SimStruct * S , int_T tid ) { real_T * lastU ;
n3qi1whofz * _rtB ; loikxjbxjg * _rtP ; ew10rzwqr2 * _rtDW ; _rtDW = ( (
ew10rzwqr2 * ) ssGetRootDWork ( S ) ) ; _rtP = ( ( loikxjbxjg * )
ssGetDefaultParam ( S ) ) ; _rtB = ( ( n3qi1whofz * ) _ssGetBlockIO ( S ) ) ;
if ( _rtDW -> nln3kpaw1s == ( rtInf ) ) { _rtDW -> nln3kpaw1s = ssGetT ( S )
; lastU = & _rtDW -> p2tj0q31lr ; } else if ( _rtDW -> ooly1cqdaf == ( rtInf
) ) { _rtDW -> ooly1cqdaf = ssGetT ( S ) ; lastU = & _rtDW -> gsdvdewgvg ; }
else if ( _rtDW -> nln3kpaw1s < _rtDW -> ooly1cqdaf ) { _rtDW -> nln3kpaw1s =
ssGetT ( S ) ; lastU = & _rtDW -> p2tj0q31lr ; } else { _rtDW -> ooly1cqdaf =
ssGetT ( S ) ; lastU = & _rtDW -> gsdvdewgvg ; } * lastU = _rtB -> pmhai3cdnu
; if ( _rtDW -> it4xhky51a == ( rtInf ) ) { _rtDW -> it4xhky51a = ssGetT ( S
) ; lastU = & _rtDW -> fcuedkogxl ; } else if ( _rtDW -> dymjgrmn00 == (
rtInf ) ) { _rtDW -> dymjgrmn00 = ssGetT ( S ) ; lastU = & _rtDW ->
dpq5ryflf5 ; } else if ( _rtDW -> it4xhky51a < _rtDW -> dymjgrmn00 ) { _rtDW
-> it4xhky51a = ssGetT ( S ) ; lastU = & _rtDW -> fcuedkogxl ; } else { _rtDW
-> dymjgrmn00 = ssGetT ( S ) ; lastU = & _rtDW -> dpq5ryflf5 ; } * lastU =
_rtB -> ohgp1nre5z ; _rtDW -> nn3yzop52g . IcNeedsLoading = 0 ; _rtDW ->
k4peu225td . IcNeedsLoading = 0 ; UNUSED_PARAMETER ( tid ) ; }
#define MDL_DERIVATIVES
static void mdlDerivatives ( SimStruct * S ) { n3qi1whofz * _rtB ; loikxjbxjg
* _rtP ; ew10rzwqr2 * _rtDW ; _rtDW = ( ( ew10rzwqr2 * ) ssGetRootDWork ( S )
) ; _rtP = ( ( loikxjbxjg * ) ssGetDefaultParam ( S ) ) ; _rtB = ( (
n3qi1whofz * ) _ssGetBlockIO ( S ) ) ; { ( ( pqmvzr1kvu * ) ssGetdX ( S ) )
-> bq3isy3pyn = _rtB -> hlobdoowm0 ; } { ( ( pqmvzr1kvu * ) ssGetdX ( S ) )
-> hnn54tmcbv = _rtB -> jcktelybde ; } { ( ( pqmvzr1kvu * ) ssGetdX ( S ) )
-> gopmpownp5 = _rtB -> mivzl3g1kg ; } { ( ( pqmvzr1kvu * ) ssGetdX ( S ) )
-> nru1t3asy3 = _rtB -> afl0c4ktu2 ; } { ( ( pqmvzr1kvu * ) ssGetdX ( S ) )
-> kfubregyho = _rtB -> cv54m4abgw ; } { ( ( pqmvzr1kvu * ) ssGetdX ( S ) )
-> nfvvmkgt4f = _rtB -> ptld5lrfqw ; } { ( ( pqmvzr1kvu * ) ssGetdX ( S ) )
-> iyebu4ash1 = _rtB -> nydpnr35ew ; } } static void mdlInitializeSizes (
SimStruct * S ) { ssSetChecksumVal ( S , 0 , 1989537487U ) ; ssSetChecksumVal
( S , 1 , 2864061156U ) ; ssSetChecksumVal ( S , 2 , 1883796930U ) ;
ssSetChecksumVal ( S , 3 , 2866438774U ) ; { mxArray * slVerStructMat = NULL
; mxArray * slStrMat = mxCreateString ( "simulink" ) ; char slVerChar [ 10 ]
; int status = mexCallMATLAB ( 1 , & slVerStructMat , 1 , & slStrMat , "ver"
) ; if ( status == 0 ) { mxArray * slVerMat = mxGetField ( slVerStructMat , 0
, "Version" ) ; if ( slVerMat == NULL ) { status = 1 ; } else { status =
mxGetString ( slVerMat , slVerChar , 10 ) ; } } mxDestroyArray ( slStrMat ) ;
mxDestroyArray ( slVerStructMat ) ; if ( ( status == 1 ) || ( strcmp (
slVerChar , "8.3" ) != 0 ) ) { return ; } } ssSetOptions ( S ,
SS_OPTION_EXCEPTION_FREE_CODE ) ; if ( ssGetSizeofDWork ( S ) != sizeof (
ew10rzwqr2 ) ) { ssSetErrorStatus ( S ,
"Unexpected error: Internal DWork sizes do "
"not match for accelerator mex file." ) ; } if ( ssGetSizeofGlobalBlockIO ( S
) != sizeof ( n3qi1whofz ) ) { ssSetErrorStatus ( S ,
"Unexpected error: Internal BlockIO sizes do "
"not match for accelerator mex file." ) ; } { int ssSizeofParams ;
ssGetSizeofParams ( S , & ssSizeofParams ) ; if ( ssSizeofParams != sizeof (
loikxjbxjg ) ) { static char msg [ 256 ] ; sprintf ( msg ,
"Unexpected error: Internal Parameters sizes do "
"not match for accelerator mex file." ) ; } } _ssSetDefaultParam ( S , (
real_T * ) & o2iu0a2jke ) ; _ssSetConstBlockIO ( S , & odcn43wyyk ) ;
rt_InitInfAndNaN ( sizeof ( real_T ) ) ; } static void
mdlInitializeSampleTimes ( SimStruct * S ) { } static void mdlTerminate (
SimStruct * S ) { }
#include "simulink.c"
