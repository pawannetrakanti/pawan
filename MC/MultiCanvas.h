#include "TString.h"
#include "TPad.h"
#include "TCanvas.h"


void makeMultiPanelCanvas(TCanvas*& canv,
                          const int columns=2,
                          const int rows=2,
                          const float leftOffset=0.5,
			  const float bottomOffset=0.5,
			  const float leftMargin=0.5,
			  const float bottomMargin=0.5,
			  const float edge=0.5) {
   if (canv==0) {
     cout<<"makeMultiPanelCanvas" << "  Got null canvas."<<endl;
     return;
   }
   canv->Clear();
   
   TPad* pad[columns][rows];

   float Xlow[columns];
   float Xup [columns];
   float Ylow[rows];
   float Yup [rows];
   float PadWidth = 
      (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
			(1.0/(1.0-edge))+(float)columns-2.0);
   float PadHeight =
      (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
			  (1.0/(1.0-edge))+(float)rows-2.0);
   Xlow[0] = leftOffset;
   Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
   Xup[columns-1] = 1;
   Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);

   Yup[0] = 1;
   Ylow[0] = 1.0-PadHeight/(1.0-edge);
   Ylow[rows-1] = bottomOffset;
   Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);

   for(int i=1;i<columns-1;i++) {
      Xlow[i] = Xup[0] + (i-1)*PadWidth;
      Xup[i] = Xup[0] + (i)*PadWidth;
   }
   int ct = 0;
   for(int i=rows-2;i>0;i--) {
      Ylow[i] = Yup[rows-1] + ct*PadHeight;
      Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
      ct++;
   }

   TString padName;
   for(int i=0;i<columns;i++) {
      for(int j=0;j<rows;j++) {
         canv->cd();
         padName = Form("p_%d_%d",i,j);
         pad[i][j] = new TPad(padName.Data(),padName.Data(),
			      Xlow[i],Ylow[j],Xup[i],Yup[j]);
         if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
         else pad[i][j]->SetLeftMargin(0);

         if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
         else pad[i][j]->SetRightMargin(0);

         if(j==0) pad[i][j]->SetTopMargin(edge);
         else pad[i][j]->SetTopMargin(0);

         if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
         else pad[i][j]->SetBottomMargin(0);

         pad[i][j]->Draw();
         pad[i][j]->cd();
         pad[i][j]->SetNumber(columns*j+i+1);
      }
   }
}
