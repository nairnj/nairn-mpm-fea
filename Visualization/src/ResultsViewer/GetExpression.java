/*
 * GetExpression
 * NairnFEAMPMViz
 * 
 * Created 1/19/2022, John A. Nairn
 * Copyright 2022, RSAC Software, All Rights Reserved
 */

import geditcom.JNFramework.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;


// The class provides a modal dialog box to run ExtractMPM
// and runs it too (all in modal window)
class GetExpression extends JNDialog
{
	private static final long serialVersionUID = 1L;
	
	private JTextField exprText=new JTextField();
	private JComboBox<PlotMenuItem> builtins=new JComboBox<PlotMenuItem>();
	DocViewer docView;
	
	//----------------------------------------------------------------------------
	// initialize
	//----------------------------------------------------------------------------

	// Create an instance of ExtractVTK dialog to pick quantites to export to VTK files
	protected GetExpression(DocViewer docCtrl,String initialExpr)
	{	// gets prompt on top and two buttons on the bottom
		super((JNFrame)docCtrl,"Expression","Enter expression to used for the next plot, or pick one from the menu"
				,"OK","Cancel",null);

		docView = docCtrl;
		
		JPanel exprPanel=new JPanel();
		GridBagLayout gridbag = new GridBagLayout();
		GridBagConstraints c = new GridBagConstraints();
		exprPanel.setLayout(gridbag);

		int top=0;
		int bottom = 2;
		int left = 2;
		int right = 2;
		c.insets=new Insets(top,left,bottom,right);			// tlbr
		c.gridx=0;
		c.gridwidth = 1;
		c.weightx = 1.0;		
		c.weighty = 1.0;
		c.fill = GridBagConstraints.HORIZONTAL;
		
		// editing field
		exprText.setText(initialExpr);
		exprText.setToolTipText("Enter any valid expression of result quantities");
		gridbag.setConstraints(exprText, c);
		exprPanel.add(exprText);
		
		// menu
		builtins.removeAllItems();
		builtins.addItem(new PlotMenuItem("Built-in expression..."));
		builtins.addItem(new PlotMenuItem("Max Princ Stress (2D):(#sxx+#syy)/2+sqrt(((#sxx-#syy)/2)^2+#sxy^2)"));
		builtins.addItem(new PlotMenuItem("Min Princ Stress (2D):(#sxx+#syy)/2-sqrt(((#sxx-#syy)/2)^2+#sxy^2)"));
		builtins.addItem(new PlotMenuItem("Max Shear Stress (2D):sqrt((#sxx-#syy)^2+2*#sxy*#sxy)"));
		builtins.addItem(new PlotMenuItem("Spherical srr (2D):#sRR*cos(#T)^2+#sZZ*sin(#T)^2+#sRZ*sin(2*#T)"));
		builtins.addItem(new PlotMenuItem("Spherical stt (2D):#sRR*sin(#T)^2+#sZZ*cos(#T)^2-#sRZ*sin(2*#T)"));
		builtins.addItem(new PlotMenuItem("Spherical szt (2D):0.5*(#sZZ-#sRR)*sin(2*#T)+#sRZ*cos(2*#T)"));
		builtins.addItem(new PlotMenuItem("Octahedral Shear (3D):sqrt((#sxx-#syy)^2+(#syy-#szz)^2+(#szz-#syy)^2+6*#sxy^2+6*#sxz^2+6*#syz^2)/3"));
		builtins.setToolTipText("Pick from any of these built-in expressions");
		builtins.setFocusable(false);
		gridbag.setConstraints(builtins, c);
		exprPanel.add(builtins);
		
		// add to dialog
		add(exprPanel,BorderLayout.CENTER);		
		
		// when use menu, insert the expression
		builtins.addItemListener(new ItemListener()
		{	public void itemStateChanged(ItemEvent e)
			{	if(builtins.getSelectedIndex()>0)
				{	PlotMenuItem pm=(PlotMenuItem)builtins.getSelectedItem();
					String value = pm.toString();
					int offset = value.indexOf(":");
					if(offset>=0)
						value = value.substring(offset+1);
					exprText.setText(value);
				}
			}
		});
		
		// initial size then override to resizable
		add(new JLabel("  "),BorderLayout.WEST);
		add(new JLabel("  "),BorderLayout.EAST);
		int width=500;
		int height=160;
		super.setSize(width,height);
		try
		{	Point windowLoc=docView.getLocationOnScreen();
			Dimension psize=docView.getSize();
			setLocation(windowLoc.x+(psize.width-width)/2,windowLoc.y+40);
		}
		catch(Exception e)
		{	setLocation(40,40);
		}
		setResizable(true);
		setVisible(true);
	}
	
	// Get the entered expression
	public String getEnteredExpression()
	{	return exprText.getText();
	}
}
		
