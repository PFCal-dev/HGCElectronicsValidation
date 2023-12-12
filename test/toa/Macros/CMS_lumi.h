#ifndef CMS_lumi_h
#define CMS_lumi_h

#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TString.h"
#include "TBox.h"
#include "TASImage.h"

//
// Global variables
//

TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold

bool writeExtraText = true; // if true, will write "Preliminary"
TString extraText   = "Preliminary";
float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.5;
float lumiTextOffset   = 0.28;
float cmsTextSize      = 0.75;
float cmsTextOffset    = 0.01;  // in case of ticks

float relPosX    = 0.045;
float relPosY    = 0.045;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.76;

bool drawLogo      = false;

void CMS_lumi( TPad* pad, int iPosX=10, TString label = "" , TString label2 = "" , bool isPreliminary = true , float lumiSize = -1 , bool printCMS = true );

#endif // #ifndef CMS_lumi_h
