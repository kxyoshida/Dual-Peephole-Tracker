import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.io.*;
import java.awt.*;
import ij.plugin.filter.*;
import ij.measure.*;
import ij.text.TextWindow;
import ij.text.TextPanel;
import java.lang.Math;



public class Dual_Peephole_Tracker implements PlugInFilter {
    ImagePlus imp;
    int w,h;
    int nSlices;
    int m=0;
    TextWindow mtw;
    TextPanel tp_opl, tp_coeff;
    int sid=0;
    int wr=9;
    int mr=4;
    int mtl=20;
    int avant=20;
    int apres=40;
    static Boolean r_to_g=true;
    static Boolean tcm=false;
    int skip;
    double cx, cy;
    double[] ax=new double[3];
    double[] bx=new double[3];
    double[] ay=new double[3];
    double[] by=new double[3];
    String title;
    String outdir = "tmp";
    /* located in the root of the ImageJ folder */

    class dPoint {
	public double x=0;
	public double y=0;
    }


	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		return DOES_16+DOES_32;
		// Updated on 25th July 2011 to allow 32-bit images
	}

	public void run(ImageProcessor ip) {
	    w=imp.getWidth()/2;
	    h=imp.getHeight();
	    nSlices=imp.getStackSize();

	    Frame[] FR = WindowManager.getNonImageWindows();

	    if (FR==null) {
		IJ.error("No tables are open.");
		return;
	    }
	    
	    int ii=0;
	    String frTitle="";
	    boolean oplWinExists=false;
	    boolean coeffWinExists=false;
	    do {
		if (FR[ii]!=null) {
		    frTitle = FR[ii].getTitle();
		    if (frTitle.toLowerCase().contains("opl")) {
			if ((FR[ii] instanceof TextWindow)) {
			    tp_opl=((TextWindow)FR[ii]).getTextPanel();
			    oplWinExists=true;
			    IJ.log("OPL file="+FR[ii].getTitle());
			} else {
			    IJ.showMessage("The opl file is not opened in the text window. Use 'import text file'.");
			    return;
			}
		    } else if (frTitle.toLowerCase().contains("coefficients")) {
			if ((FR[ii] instanceof TextWindow)) {
			    tp_coeff=((TextWindow)FR[ii]).getTextPanel();
			    coeffWinExists=true;
			    IJ.log("coefficients file="+FR[ii].getTitle());
			} else {
			    IJ.showMessage("The coefficients file is not opened in the text window. Use 'import text file'.");
			    return;
			}
		    }
		}
	    } while (++ii<FR.length);
	    
	    String line_opl=null;
	    int imax_opl=tp_opl.getLineCount();

	    int[] fr=new int[imax_opl];
	    double[] centX=new double[imax_opl];
	    double[] centY=new double[imax_opl];

	    title=imp.getShortTitle();

	    GenericDialog gd = new GenericDialog("Mask_Spot", IJ.getInstance());
	    gd.addNumericField("Window Radius:", wr, 0);
	    gd.addNumericField("Measure Window Radius:", mr, 0);
	    gd.addNumericField("Minimal Track Length", mtl, 0);
	    gd.addNumericField("Pre-track period", avant, 0);
	    gd.addNumericField("Post-track period", apres, 0);
	    gd.addCheckbox("Reverse mapping (Red to Green)", r_to_g);
	    gd.addCheckbox("Tune to the centre of mass", tcm);
	    gd.addStringField("Output directory for mini imagestacks: ", outdir);
	    gd.showDialog();
	    if (gd.wasCanceled()) 
		return;

	    wr=(int)gd.getNextNumber();
	    mr=(int)gd.getNextNumber();
	    mtl=(int)gd.getNextNumber();
	    avant=(int)gd.getNextNumber();
	    apres=(int)gd.getNextNumber();
	    r_to_g=gd.getNextBoolean();
	    tcm=gd.getNextBoolean();
	    outdir=gd.getNextString();
	    
	    IJ.log("Window Radius="+wr);
	    IJ.log("Measure Window Radius="+mr);
	    IJ.log("Minimal Track Length="+mtl);
	    IJ.log("Pre-track Period="+avant);
	    IJ.log("Post-track Period="+apres);
	    IJ.log("Reverse mapping (Red to Green)="+r_to_g);
	    IJ.log("Tune to the centre of mass="+tcm);
	    IJ.log("Output directory="+outdir);
	    
	    if (r_to_g) 
		title=title+"_RtoG";
	    else
		title=title+"_GtoR";

	    readCoefficients();
	    
	    int j=0;
	    int jmax=0;
	    
	    for (int i=0;i<imax_opl;i++) {
		line_opl=tp_opl.getLine(i);
		String[] cl=line_opl.split("\t");

		if (Integer.parseInt(cl[0])==sid) {
		    // while spot id stays the same, continue reading lines from opl file.
		    fr[j]=Integer.parseInt(cl[1]);
		    centX[j]=Double.parseDouble(cl[2]);
		    centY[j]=Double.parseDouble(cl[3]);
		    j++;
		    jmax=j;
		} else {
		    // in case spot id has changed, stop reading and go to chaseSpot routine
		    // fr, centX and centY are set for 0<=j<=jmax-1
		    if (fr[0]!=0)
			chaseSpot(ip, jmax, fr, centX, centY);
		    
		    j=0;
		    sid=Integer.parseInt(cl[0]);
		    fr[0]=Integer.parseInt(cl[1]);
		    centX[0]=Double.parseDouble(cl[2]);
		    centY[0]=Double.parseDouble(cl[3]);
		    j++;
		    jmax=j;
		}
	    }
	}


    void chaseSpot(ImageProcessor ip, int jmax, int[] fr, double[] centX, double[] centY) {
	double centXo=0;
	double centYo=0;
	int slicemin=fr[0]-avant;
	int slicemax=fr[jmax-1]+apres;
	Boolean out_of_scope=false;
	FloatProcessor mip=new FloatProcessor(4*wr+2,2*wr+1);
	ImageStack mis=new ImageStack(4*wr+2,2*wr+1);
	if (jmax>=mtl && slicemin>=1 && slicemax<=nSlices) {
	    //			IJ.log("sid,jmax,fr[0],fr[jmax-1],avant,apres="+sid+"\t"+jmax+"\t"+fr[0]+"\t"+fr[jmax-1]+"\t"+avant+"\t"+apres);
	    for (int slice=slicemin;slice<=slicemax;slice++) {
		int jj=slice-slicemin-avant;

		ip=imp.getStack().getProcessor(slice);
		
		if (slice<slicemin+avant) {
		    // for slicmin<=slice<=slicemin+avant-1 (=fr[0]-1)
		    centXo=centX[0];
		    centYo=centY[0];
		    // Do not tune centre here!
		} else if (slice>slicemax-apres) {
		    // for (fr[jmax-1]+1=) slicemax-apres+1<=slice<=slicemax
		    centXo=centX[jmax-1];
		    centYo=centY[jmax-1];
		    // Do not tune centre here!
		}  else {
		    centXo=centX[jj];
		    centYo=centY[jj];
		    int x_offset=0;
		    if (r_to_g) 
			x_offset=w;
		    else
			x_offset=0;

		    OvalRoi or=new OvalRoi(x_offset+(int)centXo-mr,(int)centYo-mr,2*mr+1,2*mr+1);
		    ip.setRoi(or);
		    ImageProcessor owip=ip.crop();
		    if (tcm) {
			dPoint o=tuneCentre(owip,centXo,centYo);
			centXo=o.x;
			centYo=o.y;
		    }
		}
		//			    IJ.log("sid,slice,jj,centXo="+sid+"\t"+slice+"\t"+jj+"\t"+centXo);
		mip=extractWindows(sid, slice,r_to_g, ip, centXo, centYo);
		if (mip!=null) 
		    mis.addSlice(title+"_sp"+sid+"_fr"+slice, mip.duplicate());
		else {
		    //				IJ.log("Out of Scope! -1");
		    out_of_scope=true;
		    break;
		}
	    }
	    if (!out_of_scope) {
		ImagePlus mimp=new ImagePlus(title+"_sp"+sid, mis);
		mimp.show();
		FileSaver fs=new FileSaver(mimp);
		fs.saveAsTiffStack(outdir+"/"+title+"_sp"+sid);
			    mimp.close();
	    } else {
		IJ.log("Out of Scope: spot id="+sid+"\tframe="+fr[0]+"-"+fr[jmax-1]+" (centXo, centYo)=("+centXo+", "+centYo+")");
	    }
	} else {
	    IJ.log("Out of Range: spot id="+sid+"\tframe="+fr[0]+"-"+fr[jmax-1]);
	}
    }

    void writeResults(int index, int frame, double rx, double ry, double gx, double gy) {
	String aLine ;
	aLine = index+"\t"+frame+"\t"+IJ.d2s(rx,4)+"\t"+IJ.d2s(ry,4)+"\t"+IJ.d2s(gx,4)+"\t"+IJ.d2s(gy,4);;
	if (mtw==null) {
	    mtw = new TextWindow(title+"_COM", "", aLine, 550, 180);
	} else
	    mtw.append(aLine);
    }

    dPoint tuneCentre(ImageProcessor owip, double centXo, double centYo) {
	ImageStatistics owis=owip.getStatistics();
	double oxd=owis.xCenterOfMass-owis.xCentroid;
	double oyd=owis.yCenterOfMass-owis.yCentroid;
	double comXo=centXo+oxd;
	double comYo=centYo+oyd;
	dPoint o=new dPoint();
	o.x=comXo;
	o.y=comYo;
	return o;
    }

    FloatProcessor extractWindows(int sid, int slice, Boolean r_to_g, ImageProcessor ip, double comXo, double comYo) {

	int x_offset;
	double comXg, comYg, comXr, comYr;
	FloatProcessor mip=new FloatProcessor(4*wr+2,2*wr+1);	
	
	double comXp=cx+ax[0]*comXo+ax[1]*Math.pow(comXo,2)+ax[2]*Math.pow(comXo,3)+bx[0]*comYo+bx[1]*Math.pow(comYo,2)+bx[2]*Math.pow(comYo,3);
	double comYp=cy+ay[0]*comXo+ay[1]*Math.pow(comXo,2)+ay[2]*Math.pow(comXo,3)+by[0]*comYo+by[1]*Math.pow(comYo,2)+by[2]*Math.pow(comYo,3);

	if (r_to_g) {
	    comXg=comXp;
	    comYg=comYp;
	    comXr=comXo;
	    comYr=comYo;
	} else {
	    comXg=comXo;
	    comYg=comYo;
	    comXr=comXp;
	    comYr=comYp;
	}

	if (comXg-wr>=0 && comXg+wr<w && comYg-wr>=0 && comYg+wr<h && comXr-wr>=0 && comXr+wr<w && comYr-wr>=0 && comYr+wr<h) {
	    for (int y=-wr;y<=wr;y++) {
		for (int x=-wr;x<=wr;x++) {
		    double fgi=ip.getBicubicInterpolatedPixel(comXg+x, comYg+y,ip);
		    double fri=ip.getBicubicInterpolatedPixel(w+comXr+x, comYr+y,ip);
		    mip.putPixelValue(wr+x,wr+y,fgi);
		    mip.putPixelValue(3*wr+1+x,wr+y,fri);
		}
	    }
	    writeResults(sid, slice, comXg, comYg, comXr, comYr);
	    return mip;
	} else
	    return null;
    }
    

    void readCoefficients() {
	    int imax_coeff=tp_coeff.getLineCount();
	    if (imax_coeff!=28) {
		IJ.showMessage("Number of coefficients does not match 2-D cubic fit.");
		return;
	    }
	    String line_coeff=null;
	    
	    if (r_to_g) 
		skip=14;
	    else
		skip=0;
	    
	    line_coeff=tp_coeff.getLine(skip);
	    cx=Double.parseDouble(line_coeff);
	    for (int i=0;i<3;i++) {
		line_coeff=tp_coeff.getLine(skip+i+1);
		ax[i]=Double.parseDouble(line_coeff);
	    }
	    for (int i=0;i<3;i++) {
		line_coeff=tp_coeff.getLine(skip+i+4);
		bx[i]=Double.parseDouble(line_coeff);
	    }
	    cy=Double.parseDouble(line_coeff);
	    for (int i=0;i<3;i++) {
		line_coeff=tp_coeff.getLine(skip+i+8);
		ay[i]=Double.parseDouble(line_coeff);
	    }
	    for (int i=0;i<3;i++) {
		line_coeff=tp_coeff.getLine(skip+i+11);
		by[i]=Double.parseDouble(line_coeff);
		//		IJ.log("by="+by[i]);
	    }
    }

}
