package FBNR;

import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.frame.*;
import java.io.*;
import javax.swing.Timer;
import java.awt.event.*;

public class progress_window extends PlugInFrame {

    private javax.swing.JProgressBar jProgressBar1;
    private javax.swing.JLabel lblMessage;
    private boolean exitStageLeft ;

	public progress_window() {
		super("Plugin_Frame");

        initComponents();

		//TextArea ta = new TextArea(15, 50);
		//add(ta);
		pack();
		GUI.center(this);
        this.setTitle("Progress");
		show();

        updateProgress();
        timer.start();
	}

    private void initComponents() {

        lblMessage = new javax.swing.JLabel();
        jProgressBar1 = new javax.swing.JProgressBar(0,100);

        lblMessage.setText("Progress...");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(
                javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(
                    layout.createParallelGroup(
                        javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jProgressBar1, 
                        javax.swing.GroupLayout.DEFAULT_SIZE,
                        277, Short.MAX_VALUE)
                    .addComponent(lblMessage)
                )
                .addContainerGap()
            )
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(
                javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(lblMessage)
                .addPreferredGap(
                    javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(jProgressBar1, 
                    javax.swing.GroupLayout.PREFERRED_SIZE, 29,
                    javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(30, Short.MAX_VALUE))
        );
    }

    void updateProgress() {
        try {
            // Open the file that is the first 
            // command line parameter
            String theFile = "" + System.getProperty("java.io.tmpdir");
            theFile += java.io.File.separator;
            theFile += "KNM_prog.txt";

            FileInputStream fstream = new FileInputStream(theFile);

            // Convert our input stream to a
            // DataInputStream
            DataInputStream in = new DataInputStream(fstream);

            // Continue to read lines while 
            // there are still some left to read
            lblMessage.setText(in.readLine());
            //set progress
            int percent = Integer.parseInt(in.readLine());
            jProgressBar1.setValue(percent);

            in.close();

            if (percent >= 100) {
                File delFile = new File(theFile);
                delFile.delete();
                exitStageLeft = closeMe();
            }
        } 
        catch (Exception e) {
            exitStageLeft = closeMe();
        }
    }

    Timer timer = new Timer(1000, new ActionListener() {
        public void actionPerformed(ActionEvent evt) {
            updateProgress();
        }
    });

    boolean closeMe(){
        IJ.wait(500);
        WindowManager.removeWindow(this);
        setVisible(false);
        timer.stop();
        dispose();

        return true;
    }

}