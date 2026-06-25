///////////////////////////////////////////////////////
//                                                   //
//  Copyright © 2005-2008 by SuperSonic Imagine, SA  //
//         Confidential - All Right Reserved         //
//                                                   //
///////////////////////////////////////////////////////

/* The Complut C function. */
/* =============================================================
 V26
 *idem V25. Correction de bug dans Reconstruct.
 V25
 * ajout du mode 0x1420504 pour accumulation de trois tirs en focalise en mode 200% 
 V24
 *modes d'accumulation 0x1410308, 0x1610308, 0x1410504, et 0x1610504 pour 3 et 4 pulses
 V21
  * Implementing photoacoustics beamforming for linear probe where transmitDistance = 0
  V20
  * Implementing
        TxAlongRx = 10 for reconstructing along the normal to line origin, and with aperture symetric along line origin
        TxAlongRx = 11 for reconstructing along the normal to line origin, and with aperture symetric along firing line origin
  V19
  * no change
  V18
  * Unchanged, but to be used with Reconstruct V18
  V17
  * Smaller gain for 12 bits ADC
  V16 2008 03 06
  * introduce modes 0210308, 0210504, 0220504
  *introduce capability to accumulate planewaves
  V15 2008 03 04
  * introduce new output mode with square root normalization for OutputTypes 100, 101, 102, and 110
  * MaxGain is reduced to keep same output intensity as for V12 with Interpolation gain set at 1024
 V14 2008 02 18
 * Correct bug when we want a force mode without immediate Lut (an immediate Lut was indeed used)
 V13 2008 02 13
 * Include new Mode management : if Mode=0, automatic mode search. If Mode >0, forced mode if coherent. If Mode=-1, accumulation mode.
 * Include 01210308, 01210504 and 01220504 modes for THI
 * Table of sound speed only with constant speed for the moment.
 * Interpolation gain set to 1024 for 10 bits ADCs
 V12 2007 09 25
 * Include new mode : 0x0010504. 5 pixels computed once in 100% mode.
 *immediate mode changed from immediate coef to immediate time but indexed coef to reduce Lut size
 *suppress unused modes
 *added possibility to input RXApod array of 4096 values when RXApodZero=0
 *store RecolutMoinsStart even when no RecoLutMonis (case of folding). This gives the size of RecolutPlus Lut
 *Add possibility to use a table of sound speeds AvSoundSpeed.
 V11 2007 04 19
 *  Change unit for PeakDelay from demod cycle to µsec.
 *  Replace structures Region and ReconInfo by AcqInfo and ReconInfo
 * Reduce max coef value from +/- 4086 to +/- 256 not to have overflow on 16 bits input.
 V10 Correct bug for TxAlongRx=0 when OutputType= 10; See around line 535
 V9 Correct for firstsample test
 V8 Correct test versus NSamples.
 V7
 *   Integration du parametre SyntAcq pour formatage des donnees en deux tirs par firing, avec donnees conjointes pour les deux tirs
 *   Gestion dun Offset memoire entre les donnees pour deux  canaux consecutifs
 *   Implementation du mode grille rectangulaire. En utilisant le code OutputType = 10 les calculs et la sortie se font selon une grille rectangulaire. Si lon est steere le PitchPix est en realite la projection verticale du pitchpix
 V6
 *   Passage de ReconDelay en µsecondes. Pour etre coherent avec TX.Delay
 V5
 *  Amelioration ponderation
    Modification de AntennaGain pour ne plus moyenner le gain dantenne par paquet de 4 pixels. Une entree de AntennaGain par pixel... Moins dartefact pour les faibles profondeurs.
    Ajout dun parametre de gain dans AntennaGain pour faire varier le gain en fonction de la profondeur. Ceci est important en flat ou la contrbution dun point est prise en compte inversement proportionnellement à la racine carrée de la profondeur. On multiplie donc en flat les pixels reconstruits proportionnellement à la racine carree de la profondeur avant de les diviser par le gain dantenne. (proportionnel a la profondeur)
    à AntennaGain a donc trois entrées par pixel : la gain dantenne pour les canaux a droite de louverture, le gain dantenne pour les canaux a gauche de louverture, et le gain global en profondeur.
 *  Troncature en mode OuputType 2
    Les shorts sont tronques en sortie à 65535
 V4  Release le 06 10 2006.
    Correction de bugs V3 :
    o   Mauvais calcul de Lut pour faibles profondeurs, lorsque tous les canaux ne contribuent pas à une region donnée.
    o   Revue de lalgorithme de ponderation en considerant separement les contributions des canaux a gauche et a droite.
    o   Revue de la taille allouée à la lut pour eviter plantage à ouverture maximale
    o   Correction de lapodisation
    Nouvelles fonctionnalites
    o   Implementation de TxAlongRx
    o   Gestion automatique debrayable du choix des modes
    o   Implementation du folding
    o   Implementation du mode 100%
    o   Implementation de lut Immediate pour petites Lut resultantes
 V3
    Cest une evolution de V2, mais avec une RecoLut qui redevient absolue, et non pas incrementale (4 fois plus grande) :

    GlobalLut est reduite aux informations de ReconMux, ReconAbsc, et ReconDelay, un entier de 4 octets  pour chaque reconstruction
    Chaque entrée de RecoLut comporte un short pour coder lapodisation, et un autre pour coder le temps.
    Louverture est symetrique par rapport à labscisse a lorigine de chaque ligne. Nouveau parametre RxAlongTx à 1
 V2
    Reduction de la taille des LUT.
    GlobalLut est reduite aux informations de ReconMux, ReconAbsc, et ReconDelay
    RecoLut code une entree en delta de temps sur un octet
    Louverture est symetrique par rapport à labscisse a lorigine de chaque ligne
    Sortie en (I,Q) float, Intensite float, ou sqrt(Intensité) en short. Nouveau parametre OutputType = 0, 1, ou 2.
 V1
    Lallocation mémoire pour le tableau de calculs intermediaires IQSum est reduite en taille dun facteur 8 (bug).
    La structure PiezoChannel est remplacée par la structure ReconMux, moins generale, mais adaptée à la structure de Mux de Rubi V1 et Rubi V2. Pas de gestion de Mux externe avec cette version.
 V0
    Algorithme de depart, livré le 15 09 2006

 * =============================================================*/

#include <math.h>

#include "complut.h"
#include <pmmintrin.h> // SSE2
#include <stdio.h>

#include "mex.h"
#include "math.h"
#include "emmintrin.h"
#include "xmmintrin.h"

//because MS compiler math.h doesn't contains the PI definition
#ifndef __GNUC__    
    #define M_PI 3.141592653589793238462643
#endif

int CompLutLinear(double *Probe, double *AcqInfo, double *ReconInfo, int *ReconMux, int *ReconAbsc, double * ReconDelay, short *Lut);
int CompLutPhased(double *Probe, double *AcqInfo, double *ReconInfo, int *ReconMux, int *ReconAbsc, double * ReconDelay, short *Lut);
int CompLutCurvedSteering(double *Probe, double *AcqInfo, double *ReconInfo, int *ReconMux, int *ReconAbsc, double * ReconDelay, short *Lut);
int CompLutLinearPhotoacoustics(double *Probe, double *AcqInfo, double *ReconInfo, int *ReconMux, int *ReconAbsc, double * ReconDelay, short *Lut);
int CompLut_FlatSynthetic_Curved(double *Probe, double *AcqInfo, double *ReconInfo, int *ReconMux, int *ReconAbsc, double * ReconDelay, short *Lut);
