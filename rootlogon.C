void loadSharedLib()
{
    // gSystem->Load("/home/john/allpix-squared/lib/libAllpixObjects.so");
    // gSystem->Load("/home/john/allpix-squared/lib/libCorryvreckanWriterObjects.so");
    gSystem->Load("lcparser/lib/liblcparser.so");
    std::cout << "Loaded lcparser library." << std::endl;
}

void globleStyle()
{
    //============================================================
    //
    //          Make graphs pretty
    //
    //============================================================
    cout << "gStyle mode requested !" << endl;

    // TGaxis::SetMaxDigits(3);

    int font = 22;

    gStyle->SetOptTitle(0);
    gStyle->SetOptDate(0);
    gStyle->SetOptStat(0);
    gStyle->SetStatColor(10);
    // gStyle->SetOptFit(0);
    gStyle->SetStatH(0.17);
    gStyle->SetStatW(0.17);
    gStyle->SetPalette(1, 0);
    gStyle->SetTextFont(font);
    gStyle->SetTextSize(0.055);
    // gStyle->SetErrorX(1);
    gStyle->SetEndErrorSize(4);
    gStyle->SetDrawBorder(0);

    gStyle->SetCanvasDefH(600);
    gStyle->SetCanvasDefW(800);
    gStyle->SetCanvasColor(10);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasBorderSize(2);
    gStyle->SetPadColor(10);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadBorderSize(0);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.10);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.02, "X");
    gStyle->SetTickLength(0.02, "Y");
    gStyle->SetPadGridX(0);
    gStyle->SetPadGridY(0);
    gStyle->SetGridColor(18);
    gStyle->SetFrameFillStyle(4000);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetFrameBorderSize(2);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(10);
    // gStyle->SetFrameLineStyle(1);

    gStyle->SetLegendFont(font);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(10);

    gStyle->SetNdivisions(510, "X");
    gStyle->SetNdivisions(510, "Y");
    gStyle->SetLabelSize(0.06, "X");
    gStyle->SetLabelSize(0.06, "Y");
    gStyle->SetLabelFont(font, "X");
    gStyle->SetLabelFont(font, "Y");
    gStyle->SetLabelOffset(0.01, "X");
    gStyle->SetLabelOffset(0.01, "Y");
    gStyle->SetTitleOffset(1.0, "X");
    gStyle->SetTitleOffset(1.0, "Y");
    gStyle->SetTitleOffset(1.0, "Z");
    gStyle->SetTitleSize(0.06, "X");
    gStyle->SetTitleSize(0.06, "Y");
    gStyle->SetTitleSize(0.06, "Z");
    gStyle->SetTitleFont(font, "X");
    gStyle->SetTitleFont(font, "Y");
    gStyle->SetTitleFont(font, "Z");
    gStyle->SetTitleColor(1);
}

// I need thick line for the frame and hist/functions
// Large white space surrounding the plot
// black font is used as the default font
void myStyle()
{
    cout << "Welcome to Shuai's style Setting !" << endl;

    TStyle *myStyle = new TStyle("myStyle", "my plots style");

    myStyle->SetPalette(1, 0);

    // use plain black on white colors
    myStyle->SetCanvasColor(10);
    myStyle->SetCanvasBorderMode(0);
    myStyle->SetCanvasBorderSize(2);
    myStyle->SetPadColor(10);
    myStyle->SetPadBorderMode(0);
    myStyle->SetPadBorderSize(0);
    myStyle->SetPadLeftMargin(0.10);
    myStyle->SetPadRightMargin(0.10);
    myStyle->SetPadTopMargin(0.08);
    myStyle->SetPadBottomMargin(0.12);
    myStyle->SetLineWidth(2); // change tick width
    myStyle->SetPadTickX(1);
    myStyle->SetPadTickY(1);
    myStyle->SetTickLength(0.02, "X");
    myStyle->SetTickLength(0.02, "Y");
    myStyle->SetPadGridX(0);
    myStyle->SetPadGridY(0);
    myStyle->SetGridColor(18);
    myStyle->SetFrameFillStyle(4000);
    myStyle->SetFrameLineWidth(2);
    myStyle->SetFrameBorderSize(2);
    myStyle->SetFrameBorderMode(0);
    myStyle->SetFrameFillColor(10);
    // gStyle->SetFrameLineStyle(1);

    // set the paper & margin sizes
    myStyle->SetPaperSize(20, 26);

    int font = 22;
    // use large Times-Roman fonts
    myStyle->SetTextFont(font);
    myStyle->SetTextSize(0.08);
    myStyle->SetLabelFont(font, "xyz");
    myStyle->SetTitleFont(font, "xyz");
    myStyle->SetLegendFont(font);
    myStyle->SetLegendBorderSize(0);
    myStyle->SetLegendFillColor(10);
    myStyle->SetStatFont(font);

    myStyle->SetLabelSize(0.04, "xyz"); // D=0.04
    myStyle->SetTitleSize(0.05, "xyz"); // D=0.02
    myStyle->SetTitleSize(0.05, "");    // main title

    myStyle->SetLabelOffset(0.015, "xyz"); // D=0.005
    myStyle->SetTitleOffset(1.0, "x");
    myStyle->SetTitleOffset(1.2, "y");
    myStyle->SetTitleOffset(1.2, "z");
    TGaxis::SetMaxDigits(3);

    // use bold lines and markers
    myStyle->SetMarkerStyle(20);
    myStyle->SetHistLineWidth(1);
    myStyle->SetLineStyleString(2, "[12 12]"); // postscript dashes

    // get rid of X error bars and y error bar caps
    // myStyle->SetErrorX(0.001);

    // do not display any of the standard histogram decorations
    myStyle->SetTitleX(0.5);
    myStyle->SetTitleAlign(23);
    // myStyle->SetTitleColor(0);
    myStyle->SetTitleStyle(0);
    myStyle->SetTitleBorderSize(0);
    myStyle->SetOptTitle(0);
    myStyle->SetOptStat(0);
    myStyle->SetOptFit(0);

    myStyle->SetStatColor(10);

    gROOT->SetStyle("myStyle");
    gROOT->ForceStyle();
}

// set color display, raibow, grayscale...
void set_color_env()
{
    cout << "Build new color environment !" << endl;

    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Int_t fcol;

    Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red[NRGBs] = {0.00, 0.00, 0.87, 1.00, 0.51};
    Double_t green[NRGBs] = {0.00, 0.81, 1.00, 0.20, 0.00};
    Double_t blue[NRGBs] = {0.51, 1.00, 0.12, 0.00, 0.00};
    fcol = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    // SetPalette has been called in the above function, color style will
    // be set according to your own rgb definition if called

    // grayscale
    /*
       double dcol = 1/double(NRGBs);
       double grey = 1;
       for(int j = 0; j < NRGBs; j++){
    // ...... Define color with RGB equal to : gray, gray, gray .......
    stops[j]=double(j)/double(NRGBs-1);
    red[j]=grey;
    blue[j]=grey;
    green[j]=grey;
    }
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    */
}

void loadFunctionMacro()
{
    // gROOT->ProcessLine(".L $HOME/Macro/function.C");
    // OR
    // gInterpreter->ProcessLine(".L $HOME/Macro/function.C");

    cout << "Self-defined functions loaded !" << endl;
}

void rootlogon()
{
    // load your system wise	root settings
    globleStyle();
    loadSharedLib();
    // myStyle();

    // loadFunctionMacro();

    // set_color_env();
}
