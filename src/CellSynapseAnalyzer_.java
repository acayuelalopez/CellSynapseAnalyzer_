import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Insets;
import java.awt.RenderingHints;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.prefs.Preferences;
import java.util.regex.PatternSyntaxException;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.DefaultComboBoxModel;
import javax.swing.DefaultListModel;
import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSlider;
import javax.swing.JSpinner;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.RowFilter;
import javax.swing.ScrollPaneConstants;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UIManager.LookAndFeelInfo;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;

import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import org.apache.commons.math3.stat.inference.TTest;
import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.plot.IntervalMarker;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.gui.ShapeRoi;
import ij.measure.Calibration;
import ij.measure.CurveFitter;
import ij.measure.Measurements;
import ij.plugin.ChannelSplitter;
import ij.plugin.CompositeConverter;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import ij.plugin.RoiRotator;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import inra.ijpb.morphology.Morphology;
import inra.ijpb.morphology.strel.SquareStrel;
import loci.plugins.in.DisplayHandler;
import loci.plugins.in.ImportProcess;
import loci.plugins.in.ImporterOptions;

public class CellSynapseAnalyzer_ implements PlugIn {

	String CELLTYPEANALYZER_IMAGES_DEFAULT_PATH;
	Preferences prefImages;
	static JLabel imageLabel;
	TextField textImages;
	JTable tableImages, tableSlice, tableS;
	DefaultTableModel modelImages, modelSlice, modelS;
	JScrollPane jScrollPaneImages, jScrollPaneSlice;
	ImagePlus[] imps, channels, arrayOfImages, impGSlices, lifs;
	ImageIcon[] icons;
	static ImagePlus imp;
	RoiManager rm;
	Roi[] roisRef, roisGTotal, roisRTotal;
	int indexMax;
	List<Double> sliceList, xrefRoiList, yrefRoiList, xnearRoiList, ynearRoiList, minDistanceList,
			minDistanceListMirror, totalRefRoisList, refArea, closestArea;
	List<Double> minDistanceListDouble;
	JCheckBox checkPhysical;
	JTextField physicalTF;
	JSpinner filterMin, filterMax;
	HistogramFilter hs2;
	IntervalMarker intervalMarker;
	ChartPanel histogram;
	JComboBox comboFilters;
	JComboBox<String> comboFitters, comboRefs, comboRefCh, comboTarCh, comboSumParam;
	DefaultListModel<String> modelList;
	JList<String> filterList;
	Double[][] dataSlice, dataS, dataFilter;
	Object[] columnHeadersS = new String[] { "N", "p-value", "Sum", "Mean", "Median", "Variance", "stdDev", "Min",
			"Max", "Q1", "Q3", "IQR" },
			columnHeadersSlice = { "Slice-ID", "Min.Distance", "xRef", "yRef", "Ref-Area", "xClosest", "yClosest",
					"Closest-Area", "Total-Rois/Slice" },
			columnNames = new Object[] { "Image", "Title", "Extension" };;
	List<Roi> roisG, roisGMirror;
	double pValue;
	double[] xG, yG, xGMirror, yGMirror;
	List<ImagePlus> impsLif = new ArrayList<ImagePlus>();
	List<ImageIcon> iconsLif = new ArrayList<ImageIcon>();
	static JTabbedPane jTabbedPane;
	KolmogorovSmirnovTest test = new KolmogorovSmirnovTest();
	List<Double> xGList, yGList;

	public void run(String arg0) {
		try {

			JFrame.setDefaultLookAndFeelDecorated(true);
			JDialog.setDefaultLookAndFeelDecorated(true);
			UIManager.setLookAndFeel("com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel");
		} catch (Exception e) {
			e.printStackTrace();
		}
		showGUI();

	}

	public void showGUI() {

		CELLTYPEANALYZER_IMAGES_DEFAULT_PATH = "images_path";
		prefImages = Preferences.userRoot();
		imageLabel = new JLabel();
		imageLabel.setBorder(BorderFactory.createTitledBorder(""));
		imageLabel.setFont(new Font(Font.DIALOG, Font.BOLD, 11));
		JButton buttonRefresh = new JButton("");
		ImageIcon iconRefresh = createImageIcon("images/refresh.png");
		Icon iconRefreshCell = new ImageIcon(iconRefresh.getImage().getScaledInstance(15, 15, Image.SCALE_SMOOTH));
		buttonRefresh.setIcon(iconRefreshCell);
		JButton buttonOpenImage = new JButton("");
		ImageIcon iconOpenImage = createImageIcon("images/openimage.png");
		Icon iconOpenImageCell = new ImageIcon(iconOpenImage.getImage().getScaledInstance(15, 15, Image.SCALE_SMOOTH));
		buttonOpenImage.setIcon(iconOpenImageCell);
		JButton buttonBrowse = new JButton("");
		ImageIcon iconBrowse = createImageIcon("images/browse.png");
		Icon iconBrowseCell = new ImageIcon(iconBrowse.getImage().getScaledInstance(15, 15, Image.SCALE_SMOOTH));
		buttonBrowse.setIcon(iconBrowseCell);
		JButton buttonProcess = new JButton("");
		ImageIcon iconProcess = createImageIcon("images/processs.png");
		Icon iconProcessCell = new ImageIcon(iconProcess.getImage().getScaledInstance(15, 15, Image.SCALE_SMOOTH));
		buttonProcess.setIcon(iconProcessCell);
		textImages = (TextField) new TextField(15);
		textImages.setText(prefImages.get(CELLTYPEANALYZER_IMAGES_DEFAULT_PATH, ""));
		DirectoryListener listenerImages = new DirectoryListener("Browse for movies...  ", textImages,
				JFileChooser.FILES_AND_DIRECTORIES);
		buttonBrowse.addActionListener(listenerImages);
		JPanel panelImagesDirect = new JPanel(new FlowLayout(FlowLayout.LEFT));
		panelImagesDirect.add(textImages);
		panelImagesDirect.add(buttonBrowse);
		JPanel panelPicture = new JPanel();
		panelPicture.setLayout(new BoxLayout(panelPicture, BoxLayout.Y_AXIS));
		JPanel bLabel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		bLabel.add(imageLabel);
		checkPhysical = new JCheckBox("Phy.Units");
		physicalTF = new JTextField(3);
		physicalTF.setEnabled(false);
		JPanel physicalPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		physicalPanel.add(checkPhysical);
		physicalPanel.add(physicalTF);
		JSeparator separator = new JSeparator(SwingConstants.VERTICAL);
		JSeparator separator1 = new JSeparator(SwingConstants.VERTICAL);
		JSeparator separator2 = new JSeparator(SwingConstants.VERTICAL);
		Dimension dime = separator.getPreferredSize();
		Dimension dime1 = separator1.getPreferredSize();
		Dimension dime2 = separator1.getPreferredSize();
		dime.height = buttonBrowse.getPreferredSize().height;
		separator.setPreferredSize(dime);
		String[] itemsRefs = new String[] { "All Refs", "&& Refs" };
		comboRefs = new JComboBox<String>();
		for (int i = 0; i < itemsRefs.length; i++)
			comboRefs.addItem(itemsRefs[i]);
		comboRefs.setSelectedIndex(1);
		comboRefs.setPreferredSize(new Dimension(80, 25));
		JLabel refChLabel = new JLabel("Reference-Ch: ");
		JLabel targChLabel = new JLabel("Target-Ch: ");
		comboRefCh = new JComboBox<String>();
		comboRefCh.setPreferredSize(new Dimension(80, 25));
		comboTarCh = new JComboBox<String>();
		comboTarCh.setPreferredSize(new Dimension(80, 25));
		JPanel refPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		refPanel.add(refChLabel);
		refPanel.add(comboRefCh);
		JPanel tarPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		tarPanel.add(targChLabel);
		tarPanel.add(comboTarCh);
		JPanel boxPanel = new JPanel();
		boxPanel.setLayout(new BoxLayout(boxPanel, BoxLayout.Y_AXIS));
		boxPanel.add(refPanel);
		boxPanel.add(tarPanel);
		dime1.height = boxPanel.getPreferredSize().height;
		dime2.height = boxPanel.getPreferredSize().height;
		separator1.setPreferredSize(dime1);
		separator2.setPreferredSize(dime2);
		bLabel.add(separator);
		bLabel.add(panelImagesDirect);
		bLabel.add(buttonRefresh);
		bLabel.add(buttonOpenImage);
		JPanel panelSB = new JPanel();
		panelSB.add(CellSynapseProcessing_.pb);
		bLabel.add(panelSB);
		tableImages = new JTable();
		modelImages = new DefaultTableModel();
		modelImages.setColumnIdentifiers(columnNames);
		tableImages.setModel(modelImages);
		jScrollPaneImages = new JScrollPane(tableImages);
		jScrollPaneImages.setPreferredSize(new Dimension(650, 200));
		modelImages = new DefaultTableModel(columnNames, 1) {

			@Override
			public Class<?> getColumnClass(int column) {
				if (getRowCount() > 0) {
					Object value = getValueAt(0, column);
					if (value != null) {
						return getValueAt(0, column).getClass();
					}
				}

				return super.getColumnClass(column);
			}

			public boolean isCellEditable(int row, int col) {
				return false;
			}

		};
		tableImages.setModel(modelImages);
		tableImages.setSelectionBackground(new Color(229, 255, 204));
		tableImages.setSelectionForeground(new Color(0, 102, 0));
		DefaultTableCellRenderer centerRenderer = new DefaultTableCellRenderer();
		centerRenderer.setHorizontalAlignment(JLabel.CENTER);
		tableImages.setDefaultRenderer(String.class, centerRenderer);
		tableImages.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		tableImages.setRowHeight(60);
		tableImages.setAutoCreateRowSorter(true);
		tableImages.getTableHeader().setDefaultRenderer(new SimpleHeaderRenderer());
		tableImages.getColumnModel().getColumn(0).setPreferredWidth(100);
		tableImages.getColumnModel().getColumn(1).setPreferredWidth(450);
		tableImages.getColumnModel().getColumn(2).setPreferredWidth(100);
		tableSlice = new JTable();
		modelSlice = new DefaultTableModel();
		modelSlice.setColumnIdentifiers(columnHeadersSlice);
		tableSlice.setModel(modelSlice);
		jScrollPaneSlice = new JScrollPane(tableSlice);
		jScrollPaneSlice.setPreferredSize(new Dimension(650, 170));
		modelSlice = new DefaultTableModel(columnHeadersSlice, 1) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public Class<?> getColumnClass(int column) {
				if (getRowCount() > 0) {
					Object value = getValueAt(0, column);
					if (value != null) {
						return getValueAt(0, column).getClass();
					}
				}

				return super.getColumnClass(column);
			}

			public boolean isCellEditable(int row, int col) {
				return false;
			}

		};
		tableSlice.setModel(modelSlice);
		tableSlice.setSelectionBackground(new Color(229, 255, 204));
		tableSlice.setSelectionForeground(new Color(0, 102, 0));
		// tableIndividual.setDefaultRenderer(ImageIcon.class, centerRenderer);
		tableSlice.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		tableSlice.setRowHeight(60);
		tableSlice.setAutoCreateRowSorter(true);
		tableSlice.getTableHeader().setDefaultRenderer(new SimpleHeaderRenderer());
		for (int u = 0; u < tableSlice.getColumnCount(); u++)
			tableSlice.getColumnModel().getColumn(u).setPreferredWidth(150);
		tableS = new JTable();
		modelS = new DefaultTableModel();
		modelS.setColumnIdentifiers(columnHeadersS);
		tableS.setModel(modelS);
		modelS = new DefaultTableModel(columnHeadersS, 1) {

			@Override
			public Class<?> getColumnClass(int column) {
				if (getRowCount() != 0) {
					Object value = getValueAt(0, column);
					if (value != null) {
						return getValueAt(0, column).getClass();
					}
				}
				return super.getColumnClass(column);
			}

			public boolean isCellEditable(int row, int col) {
				return false;
			}

		};
		tableS.setModel(modelS);
		tableS.setSelectionBackground(new Color(229, 255, 204));
		tableS.setSelectionForeground(new Color(0, 102, 0));
		tableS.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		tableS.setRowHeight(60);
		tableS.setAutoCreateRowSorter(true);
		tableS.getTableHeader().setDefaultRenderer(new SimpleHeaderRenderer());
		for (int u = 0; u < tableS.getColumnCount(); u++)
			tableS.getColumnModel().getColumn(u).setPreferredWidth(150);
		JScrollPane jScrollPaneS = new JScrollPane(tableS);
		jScrollPaneS.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		jScrollPaneS.setPreferredSize(new Dimension(650, 83));
		JPanel mainPanel = new JPanel();
		mainPanel.setLayout(new BoxLayout(mainPanel, BoxLayout.Y_AXIS));
		JPanel imagePanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		imagePanel.add(jScrollPaneImages);
		JPanel slicePanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		slicePanel.add(jScrollPaneSlice);
		JPanel sPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		sPanel.add(jScrollPaneS);
		panelPicture.add(Box.createVerticalStrut(25));
		panelPicture.add(slicePanel);
		panelPicture.add(Box.createVerticalStrut(3));
		JLabel sumParameterLabel = new JLabel("Summary-Parameter: ");
		comboSumParam = new JComboBox<String>();
		JPanel sumParamPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		sumParamPanel.add(sumParameterLabel);
		sumParamPanel.add(comboSumParam);
		JPanel lastPanel0 = new JPanel(new FlowLayout(FlowLayout.LEFT));
		lastPanel0.add(comboRefs);
		lastPanel0.add(physicalPanel);
		lastPanel0.add(buttonProcess);
		JPanel lastPanel1 = new JPanel();
		lastPanel1.setLayout(new BoxLayout(lastPanel1, BoxLayout.Y_AXIS));
		lastPanel1.add(lastPanel0);
		lastPanel1.add(sumParamPanel);
		JPanel lastPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		lastPanel.add(separator1);
		lastPanel.add(boxPanel);
		lastPanel.add(separator2);
		lastPanel.add(lastPanel1);
		panelPicture.add(lastPanel);
		panelPicture.add(sPanel);
		mainPanel.add(imagePanel);
		mainPanel.add(bLabel);
		JPanel filtersMin = new JPanel(new FlowLayout(FlowLayout.CENTER));
		filterMin = new JSpinner(new SpinnerNumberModel(30, 0, 5000, 1));
		filterMin.setPreferredSize(new Dimension(60, 20));
		JSlider sliderMin = new JSlider(0, 300, 50);
		sliderMin.setPreferredSize(new Dimension(150, 15));
		JLabel filterMinLabel = new JLabel(" Min :  ");
		filtersMin.add(filterMinLabel);
		filtersMin.add(sliderMin);
		filtersMin.add(Box.createHorizontalStrut(2));
		filtersMin.add(filterMin);

		JPanel filtersMax = new JPanel(new FlowLayout(FlowLayout.CENTER));
		filterMax = new JSpinner(new SpinnerNumberModel(200, 0, 5000, 1));
		filterMax.setPreferredSize(new Dimension(60, 20));
		JSlider sliderMax = new JSlider(0, 300, 150);
		sliderMax.setPreferredSize(new Dimension(150, 15));
		JLabel filterMaxLabel = new JLabel(" Max :  ");
		filtersMax.add(filterMaxLabel);
		filtersMax.add(sliderMax);
		filtersMax.add(Box.createHorizontalStrut(2));
		filtersMax.add(filterMax);

		JPanel boxPanel2 = new JPanel();
		boxPanel2.setLayout(new BoxLayout(boxPanel2, BoxLayout.Y_AXIS));
		JPanel chartPanel2 = new JPanel();
		hs2 = new HistogramFilter();
		intervalMarker = new IntervalMarker(0, 0, new Color(229, 255, 204), new BasicStroke(), new Color(0, 102, 0),
				new BasicStroke(1.5f), 0.5f);
		histogram = hs2.createChartPanel("Distance", new double[] { 0.0, 0.0, 0.0 }, 100, intervalMarker);
		chartPanel2.add(histogram);
		boxPanel2.add(chartPanel2);
		JPanel controlPanel2 = hs2.createControlPanel();
		boxPanel2.add(controlPanel2);

		JPanel filtersMain = new JPanel();
		filtersMain.setLayout(new BoxLayout(filtersMain, BoxLayout.Y_AXIS));
		filtersMain.add(boxPanel2);
		filtersMain.add(filtersMin);
		filtersMain.add(filtersMax);

		String itemFilters[] = { "Slice-ID", "Min.Distance" };

		comboFilters = new JComboBox<String>();
		for (int i = 0; i < columnHeadersSlice.length; i++)
			comboFilters.addItem(columnHeadersSlice[i]);
		comboFilters.setSelectedIndex(0);
		comboFilters.setOpaque(true);
		comboFilters.setPreferredSize(new Dimension(80, 25));
		JButton btnAdd = new JButton();
		ImageIcon iconAdd = createImageIcon("images/plus.png");
		Icon addCell = new ImageIcon(iconAdd.getImage().getScaledInstance(15, 15, Image.SCALE_SMOOTH));
		btnAdd.setIcon(addCell);
		btnAdd.setToolTipText("Click this button to add a filter to list.");
		JButton btnRem = new JButton();
		ImageIcon iconRem = createImageIcon("images/minus.png");
		Icon remCell = new ImageIcon(iconRem.getImage().getScaledInstance(15, 15, Image.SCALE_SMOOTH));
		btnRem.setIcon(remCell);
		btnRem.setToolTipText("Click this button to remove a filter from list.");
		JButton btnCF = new JButton();
		ImageIcon iconCF = createImageIcon("images/curvefitter.png");
		Icon cfCell = new ImageIcon(iconCF.getImage().getScaledInstance(15, 15, Image.SCALE_SMOOTH));
		btnCF.setIcon(cfCell);
		btnCF.setToolTipText("Click this button to get distribution.");
		modelList = new DefaultListModel<String>();
		filterList = new JList<String>(modelList);
		JPanel boxMPButtons = new JPanel();
		boxMPButtons.setLayout(new BoxLayout(boxMPButtons, BoxLayout.Y_AXIS));
		boxMPButtons.add(btnAdd);
		boxMPButtons.add(btnRem);

		JPanel listButtons = new JPanel(new FlowLayout(FlowLayout.LEFT));
		JScrollPane scrollFilterList = new JScrollPane(filterList);
		scrollFilterList.setPreferredSize(new Dimension(150, 100));
		listButtons.add(scrollFilterList);
		listButtons.add(boxMPButtons);

		JPanel filterFeature = new JPanel(new FlowLayout(FlowLayout.LEFT));
		filterFeature.add(comboFilters);
		JPanel boxFilter = new JPanel();
		boxFilter.setLayout(new BoxLayout(boxFilter, BoxLayout.Y_AXIS));
		JButton csvButton = new JButton("");
		ImageIcon iconCsv = createImageIcon("images/csv.png");
		Icon csvCell = new ImageIcon(iconCsv.getImage().getScaledInstance(15, 15, Image.SCALE_SMOOTH));
		csvButton.setIcon(csvCell);
		csvButton.setToolTipText("Click this button to export table data as csv.file.");
		JButton filterButton = new JButton("");
		ImageIcon iconFilter = createImageIcon("images/filter.png");
		Icon filterCell = new ImageIcon(iconFilter.getImage().getScaledInstance(15, 15, Image.SCALE_SMOOTH));
		filterButton.setIcon(filterCell);
		filterButton.setToolTipText("Click this button to filter your spot detections.");
		JButton resetFilterButton = new JButton("");
		ImageIcon iconResetFilter = createImageIcon("images/resetfilter.png");
		Icon resetFilterCell = new ImageIcon(iconResetFilter.getImage().getScaledInstance(15, 15, Image.SCALE_SMOOTH));
		resetFilterButton.setIcon(resetFilterCell);
		resetFilterButton.setToolTipText("Click this button to reset your filtering selection.");
		boxFilter.add(filterFeature);
		boxFilter.add(Box.createVerticalStrut(1));
		boxFilter.add(listButtons);
		boxFilter.add(Box.createVerticalStrut(3));
		JPanel buttonFilter = new JPanel(new FlowLayout(FlowLayout.CENTER));
		buttonFilter.add(filterButton);
		buttonFilter.add(resetFilterButton);
		boxFilter.add(buttonFilter);
		JSeparator separator3 = new JSeparator(SwingConstants.HORIZONTAL);
		Dimension dime3 = separator3.getPreferredSize();
		dime3.width = boxFilter.getPreferredSize().width;
		separator3.setPreferredSize(dime3);
		boxFilter.add(Box.createVerticalStrut(3));
		boxFilter.add(separator3);
		JPanel panelCF = new JPanel(new FlowLayout(FlowLayout.LEFT));
		panelCF.add(btnCF);
		String[] itemsFitters = new String[] { "Linear", "2Deg.Poly", "3Deg.Poly", "24Deg.Poly", "5Deg.Poly",
				"6Deg.Poly", "7Deg.Poly", "8Deg.Poly", "Power", "Exponential", "Log", "Gaussian" };
		comboFitters = new JComboBox<String>();
		for (int i = 0; i < itemsFitters.length; i++)
			comboFitters.addItem(itemsFitters[i]);
		Dimension size = new Dimension();
		size.width = comboFilters.getPreferredSize().width + 5;
		panelCF.add(comboFitters);
		boxFilter.add(Box.createVerticalStrut(3));
		boxFilter.add(panelCF);
		JPanel csvPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		csvPanel.add(csvButton);
		boxFilter.add(csvPanel);
		JPanel filterDef = new JPanel(new FlowLayout(FlowLayout.LEFT));
		filterDef.add(filtersMain);
		filterDef.add(boxFilter);
		panelPicture.add(filterDef);
		jTabbedPane = new JTabbedPane();
		CellSynapseProcessing_.pb.setString("...");
		CellSynapseProcessing_ csA = new CellSynapseProcessing_();
		csA.showGUI();
		mainPanel.add(jTabbedPane);

		JFrame frame = new JFrame();
		frame.setTitle("Cell-Synapse Analyzer");
		frame.setResizable(true);
		frame.add(mainPanel);
		frame.pack();
		frame.setSize(660, 980);
		frame.setLocationRelativeTo(null);
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.setVisible(true);

		sliderMin.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {

				filterMin.setValue(sliderMin.getValue());
				intervalMarker.setStartValue(sliderMin.getValue());

			}
		});

		filterMin.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
				sliderMin.setValue((int) filterMin.getValue());
				intervalMarker.setStartValue((int) filterMin.getValue());

			}
		});

		sliderMax.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
				filterMax.setValue(sliderMax.getValue());
				intervalMarker.setEndValue(sliderMax.getValue());
			}
		});

		filterMax.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent e) {
				sliderMax.setValue((int) filterMax.getValue());
				intervalMarker.setEndValue((int) filterMax.getValue());

			}
		});
		filterButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {

				filterTool();

			}
		});
		resetFilterButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {

				resetFilter();

			}
		});
		csvButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {

				csvExport();

			}
		});

		btnAdd.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {

				List<String> listFilters = new ArrayList<String>();

				if (filterList.getModel().getSize() < 1)
					modelList.addElement((String) comboFilters.getSelectedItem() + ":  (" + filterMin.getValue() + ","
							+ filterMax.getValue() + ")");

				if (filterList.getModel().getSize() >= 1) {
					for (int i = 0; i < filterList.getModel().getSize(); i++)
						listFilters.add(String.valueOf(filterList.getModel().getElementAt(i).substring(0,
								filterList.getModel().getElementAt(i).lastIndexOf(":"))));

					if (listFilters.contains(comboFilters.getSelectedItem().toString()) == false)
						modelList.addElement((String) comboFilters.getSelectedItem() + ":  (" + filterMin.getValue()
								+ "," + filterMax.getValue() + ")");

					if (listFilters.contains(comboFilters.getSelectedItem().toString()) == true)
						return;

				}

			}
		});
		btnRem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {

				try {
					int[] indexes = filterList.getSelectedIndices();
					for (int i = 0; i < indexes.length; i++)
						modelList.remove(indexes[i]);
				} catch (Exception e1) {
					e1.printStackTrace();
				}

			}
		});

		checkPhysical.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent e) {
				if (e.getStateChange() == ItemEvent.SELECTED)
					physicalTF.setEnabled(true);
				if (e.getStateChange() == ItemEvent.DESELECTED)
					physicalTF.setEnabled(false);
			}
		});
		comboFilters.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				String selectedName = (String) comboFilters.getSelectedItem();
				int selectedIndex = (int) comboFilters.getSelectedIndex();

				double values[] = new double[dataSlice.length];
				for (int r = 0; r < dataSlice.length; r++)
					values[r] = dataSlice[r][selectedIndex];

				int i;
				double maxValue = values[0];
				double minValue = values[0];
				for (i = 1; i < values.length; i++) {
					if (values[i] > maxValue)
						maxValue = values[i];
					if (values[i] < minValue)
						minValue = values[i];
				}
				sliderMin.setMinimum((int) minValue - 1);
				sliderMin.setMaximum((int) maxValue + 1);
				sliderMax.setMinimum((int) minValue - 1);
				sliderMax.setMaximum((int) maxValue + 1);

				hs2.addHistogramSeries(selectedName, values, dataSlice.length, intervalMarker);

			}
		});
		buttonRefresh.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				refreshAction();
			}
		});

		btnCF.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				curveFitter();
			}
		});

		buttonOpenImage.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				List<ImagePlus> impAnalClose = new ArrayList<ImagePlus>();
				int[] IDs = WindowManager.getIDList();
				if (IDs != null)
					for (int i = 0; i < IDs.length; i++)
						impAnalClose.add(WindowManager.getImage(IDs[i]));

				if (tableImages.getSelectedRow() != -1) {
					if (IDs != null)
						for (int i = 0; i < IDs.length; i++)
							impAnalClose.get(i).hide();
					if (impsLif.isEmpty() == false)
						imp = impsLif.get(tableImages.getSelectedRow());
					if (impsLif.isEmpty() == true)
						imp = imps[tableImages.getSelectedRow()];
				}
				if (imp == null)
					IJ.error("Please, select an image within the main directory.");
				if (imp != null) {
					if (imp.getTitle().contains(".lif") == true) {
						imp = CompositeConverter.makeComposite(imp).duplicate();
						imp.setDisplayMode(IJ.COLOR);
					}
					if (imp.getTitle().contains(".lif") == false)
						imp.setDisplayMode(IJ.COLOR);
					imp.show();
					String title = imp.getTitle().substring(0, imp.getTitle().lastIndexOf("."));
					imageLabel.setText(title);

					comboRefCh.removeAllItems();
					comboTarCh.removeAllItems();
					if (imp.getTitle().contains(".lif") == Boolean.TRUE) {
						for (int i = 0; i < imp.getChannelProcessor().getNChannels(); i++) {
							comboRefCh.addItem("Ch-" + (i + 1));
							comboTarCh.addItem("Ch-" + (+i + 1));
						}
						CellSynapseProcessing_.comboCellCh.removeAllItems();
						CellSynapseProcessing_.comboSynCh.removeAllItems();
						CellSynapseProcessing_.comboProtCh.removeAllItems();
						for (int i = 0; i < imp.getNChannels(); i++) {
							CellSynapseProcessing_.comboCellCh.addItem("Ch-" + (i + 1));
							CellSynapseProcessing_.comboSynCh.addItem("Ch-" + (+i + 1));
							CellSynapseProcessing_.comboProtCh.addItem("Ch-" + (+i + 1));
						}
					}
					if (imp.getTitle().contains(".lsm") == Boolean.TRUE
							|| imp.getTitle().contains(".czi") == Boolean.TRUE) {
						for (int i = 0; i < imp.getChannelProcessor().getNChannels(); i++) {
							comboRefCh.addItem("Ch-" + (i + 1));
							comboTarCh.addItem("Ch-" + (+i + 1));
						}
						CellSynapseProcessing_.comboCellCh.removeAllItems();
						CellSynapseProcessing_.comboSynCh.removeAllItems();
						CellSynapseProcessing_.comboProtCh.removeAllItems();
						for (int i = 0; i < imp.getNChannels(); i++) {
							CellSynapseProcessing_.comboCellCh.addItem("Ch-" + (i + 1));
							CellSynapseProcessing_.comboSynCh.addItem("Ch-" + (+i + 1));
							CellSynapseProcessing_.comboProtCh.addItem("Ch-" + (+i + 1));
						}
					}

				}
			}
		});
		tableSlice.addMouseListener(new MouseAdapter() {

			@Override
			public void mouseReleased(final MouseEvent e) {
				int indexSelected = tableSlice.getSelectedRow();
				rm.select(indexSelected);

			}
		});
	}

	public static Image getScaledImage(Image srcImg, int w, int h) {
		BufferedImage resizedImg = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
		Graphics2D g2 = resizedImg.createGraphics();
		g2.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
		g2.drawImage(srcImg, 0, 0, w, h, null);
		g2.dispose();
		return resizedImg;
	}

	public void csvExport() {
		JFrame parentFrame = new JFrame();

		JFileChooser fileChooser = new JFileChooser();
		fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		fileChooser.setDialogTitle("Select a Directory to save data");
		int userSelection = fileChooser.showSaveDialog(parentFrame);
		File fileToSave = null;
		if (userSelection == JFileChooser.APPROVE_OPTION) {
			fileToSave = fileChooser.getSelectedFile();

			try {

				TableModel modelSlice = tableSlice.getModel();
				TableModel modelSum = tableS.getModel();
				FileWriter csv = new FileWriter(
						new File(fileToSave.getAbsolutePath() + File.separator + imageLabel.getText() + "_data.csv"));

				for (int i = 0; i < modelSlice.getColumnCount(); i++) {
					csv.write(modelSlice.getColumnName(i) + ",");
				}

				csv.write("\n");

				for (int i = 0; i < tableSlice.getRowCount(); i++) {
					for (int j = 0; j < tableSlice.getColumnCount(); j++) {

						csv.write(tableSlice.getModel().getValueAt(tableSlice.convertRowIndexToModel(i), j).toString()
								+ ",");
					}
					csv.write("\n");
				}
				csv.write("\n");
				csv.write("\n");
				csv.write("\n");
				for (int i = 0; i < modelSum.getColumnCount(); i++) {
					csv.write(modelSum.getColumnName(i) + ",");
				}

				csv.write("\n");

				for (int i = 0; i < modelSum.getRowCount(); i++) {
					for (int j = 0; j < modelSum.getColumnCount(); j++) {
						csv.write(modelSum.getValueAt(i, j).toString() + ",");
					}
					csv.write("\n");

				}

				csv.close();

			} catch (IOException e) {
				e.printStackTrace();
			}
		}

	}

	public void curveFitter() {

		double[] xValues = null;
		double[] yValues = null;
		List<Double> xValuesL = new ArrayList<Double>();
		List<Double> yValuesL = new ArrayList<Double>();

		if (tableSlice.getRowCount() != dataSlice.length) {
			for (int row = 0; row < tableSlice.getRowCount(); row++) {
				xValuesL.add((double) row + 1);
				yValuesL.add(((Double) tableSlice.getModel().getValueAt(tableSlice.convertRowIndexToModel(row), 1)));
			}
			xValues = new double[xValuesL.size()];
			yValues = new double[yValuesL.size()];
			for (int i = 0; i < xValuesL.size(); i++) {
				xValues[i] = xValuesL.get(i).doubleValue();
				yValues[i] = yValuesL.get(i).doubleValue();
			}

		}
		if (tableSlice.getRowCount() == dataSlice.length) {
			for (int row = 0; row < modelSlice.getRowCount(); row++) {
				xValuesL.add((double) row + 1);
				yValuesL.add(((Double) modelSlice.getValueAt(tableSlice.convertRowIndexToModel(row), 1)));

			}
			xValues = new double[xValuesL.size()];
			yValues = new double[yValuesL.size()];
			for (int i = 0; i < xValuesL.size(); i++) {
				xValues[i] = xValuesL.get(i).doubleValue();
				yValues[i] = yValuesL.get(i).doubleValue();
			}
		}
		// Plot plot = new Plot("Curve Fitter", "x", "y");
		// plot.add("circle", xValues, yValues);

		CurveFitter cf = new CurveFitter(xValues, yValues);

		if (comboFitters.getSelectedItem().equals("Linear"))
			cf.doFit(CurveFitter.STRAIGHT_LINE);
		if (comboFitters.getSelectedItem().equals("2Deg.Poly"))
			cf.doFit(CurveFitter.POLY2);
		if (comboFitters.getSelectedItem().equals("3Deg.Poly"))
			cf.doFit(CurveFitter.POLY3);
		if (comboFitters.getSelectedItem().equals("4Deg.Poly"))
			cf.doFit(CurveFitter.POLY4);
		if (comboFitters.getSelectedItem().equals("5Deg.Poly"))
			cf.doFit(CurveFitter.POLY5);
		if (comboFitters.getSelectedItem().equals("6Deg.Poly"))
			cf.doFit(CurveFitter.POLY6);
		if (comboFitters.getSelectedItem().equals("7Deg.Poly"))
			cf.doFit(CurveFitter.POLY7);
		if (comboFitters.getSelectedItem().equals("8Deg.Poly"))
			cf.doFit(CurveFitter.POLY8);
		if (comboFitters.getSelectedItem().equals("Power"))
			cf.doFit(CurveFitter.POWER);
		if (comboFitters.getSelectedItem().equals("Exponential"))
			cf.doFit(CurveFitter.EXPONENTIAL);
		if (comboFitters.getSelectedItem().equals("Log"))
			cf.doFit(CurveFitter.LOG);
		if (comboFitters.getSelectedItem().equals("Gaussian"))
			cf.doFit(CurveFitter.GAUSSIAN);
		Plot plot = cf.getPlot();
		plot.addLabel(0, 0.35, "sdtDev = " + IJ.d2s(cf.getSD(), 4));
		plot.setXYLabels("Roi", "Min.Distance");
		plot.show();

	}

	public void filterTool() {
		if (filterList.getModel().getSize() == 0) {
			IJ.error("You should add at least one parameter for filtering.");
			return;
		}
		if (filterList.getModel().getSize() >= 1) {
			List<String> listFiltersName = new ArrayList<String>();
			List<String> listFiltersMin = new ArrayList<String>();
			List<String> listFiltersMax = new ArrayList<String>();
			for (int i = 0; i < filterList.getModel().getSize(); i++) {
				listFiltersName.add(String.valueOf(filterList.getModel().getElementAt(i).substring(0,
						filterList.getModel().getElementAt(i).lastIndexOf(":"))));
				listFiltersMin.add(String.valueOf(filterList.getModel().getElementAt(i).substring(
						filterList.getModel().getElementAt(i).lastIndexOf("("),
						filterList.getModel().getElementAt(i).lastIndexOf(","))).replace("(", ""));
				listFiltersMax.add(String.valueOf(filterList.getModel().getElementAt(i).substring(
						filterList.getModel().getElementAt(i).lastIndexOf(","),
						filterList.getModel().getElementAt(i).lastIndexOf(")"))).replace(",", ""));

			}

			DefaultComboBoxModel<String> modelComboFilters = (DefaultComboBoxModel<String>) comboFilters.getModel();
			List<Integer> indexesComboFilters = new ArrayList<Integer>();
			for (int i = 0; i < listFiltersName.size(); i++)
				indexesComboFilters.add(modelComboFilters.getIndexOf(listFiltersName.get(i)));
			TableRowSorter rowSorterR = new TableRowSorter<>(tableSlice.getModel());
			tableSlice.setRowSorter(rowSorterR);
			List<RowFilter<TableModel, Integer>> low = new ArrayList<>();
			List<RowFilter<TableModel, Integer>> high = new ArrayList<>();
			for (int i = 0; i < indexesComboFilters.size(); i++) {

				low.add(RowFilter.numberFilter(RowFilter.ComparisonType.AFTER,
						Double.parseDouble(listFiltersMin.get(i)),
						tableSlice.convertColumnIndexToModel(indexesComboFilters.get(i))));
				high.add(RowFilter.numberFilter(RowFilter.ComparisonType.BEFORE,
						Double.parseDouble(listFiltersMax.get(i)),
						tableSlice.convertColumnIndexToModel(indexesComboFilters.get(i))));
			}

			List<RowFilter<TableModel, Integer>> listOfFilters = new ArrayList<>();
			for (int i = 0; i < low.size(); i++)
				listOfFilters.add(low.get(i));
			for (int i = 0; i < high.size(); i++)
				listOfFilters.add(high.get(i));

			final RowFilter<TableModel, Integer> filter = RowFilter.andFilter(listOfFilters);

			for (int i = 0; i < indexesComboFilters.size(); i++)
				if (Double.parseDouble(listFiltersMin.get(i)) == 0 && Double.parseDouble(listFiltersMax.get(i)) == 0) {
					rowSorterR.setRowFilter(null);
				} else {
					try {
						rowSorterR.setRowFilter(filter);
					} catch (PatternSyntaxException pse) {
						System.out.println("Bad regex pattern");
					}
				}

		}

		List<Integer> indexFilter = new ArrayList<Integer>();
		List<Double> minDistFil = new ArrayList<Double>();
		List<Double> minDistMirrorFil = new ArrayList<Double>();
		List<Double> listFiltered = new ArrayList<Double>();

		for (int row = 0; row < tableSlice.getRowCount(); row++) {
			indexFilter.add(tableSlice.convertRowIndexToModel(row));
			listFiltered.add((Double) tableSlice.getModel().getValueAt(tableSlice.convertRowIndexToModel(row),
					comboSumParam.getSelectedIndex() + 1));
		}

		for (int row = 0; row < indexFilter.size(); row++) {
			minDistFil.add(minDistanceList.get(indexFilter.get(row)));
			minDistMirrorFil.add(minDistanceListMirror.get(indexFilter.get(row)));

		}
		double[] minDistanceArrayFil = new double[minDistFil.size()];
		double[] minDistanceMirrorArrayFil = new double[minDistMirrorFil.size()];

		for (int x = 0; x < minDistanceArrayFil.length; x++) {
			minDistanceArrayFil[x] = minDistFil.get(x).doubleValue();
			minDistanceMirrorArrayFil[x] = minDistMirrorFil.get(x).doubleValue();

		}

		double pValueFil = test.kolmogorovSmirnovTest(minDistanceArrayFil, minDistanceMirrorArrayFil, true);
		modelS.setValueAt((double) (double) tableSlice.getRowSorter().getViewRowCount(),
				tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(0));
		modelS.setValueAt(pValueFil, tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(1));
		modelS.setValueAt(Math.round(listFiltered.stream().mapToDouble(Double::doubleValue).sum() * 1000.0) / 1000.0,
				tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(2));
		modelS.setValueAt(
				Math.round(listFiltered.stream().mapToDouble(Double::doubleValue).sum() / listFiltered.size() * 1000.0)
						/ 1000.0,
				tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(3));
		modelS.setValueAt(Math.round(CellSynapseProcessing_.median(listFiltered) * 1000.0) / 1000.0,
				tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(4));
		modelS.setValueAt(Math.round(CellSynapseProcessing_.variance(listFiltered) * 1000.0) / 1000.0,
				tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(5));
		modelS.setValueAt(Math.round(CellSynapseProcessing_.sd(listFiltered) * 1000.0) / 1000.0,
				tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(6));
		modelS.setValueAt(Math.round(Collections.min(listFiltered) * 1000.0) / 1000.0, tableS.convertRowIndexToModel(0),
				tableS.convertColumnIndexToModel(7));
		modelS.setValueAt(Math.round(Collections.max(listFiltered) * 1000.0) / 1000.0, tableS.convertRowIndexToModel(0),
				tableS.convertColumnIndexToModel(8));
		modelS.setValueAt(Math.round(CellSynapseProcessing_.Q1(listFiltered, listFiltered.size()) * 1000.0) / 1000.0,
				tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(9));
		modelS.setValueAt(Math.round(CellSynapseProcessing_.Q3(listFiltered, listFiltered.size()) * 1000.0) / 1000.0,
				tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(10));
		modelS.setValueAt(Math.round(CellSynapseProcessing_.IQR(listFiltered, listFiltered.size()) * 1000.0) / 1000.0,
				tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(11));

		for (int j = 0; j < indexFilter.size(); j++) {
			roisG.get(indexFilter.get(j).intValue())
					.setPosition(sliceList.get(indexFilter.get(j).intValue()).intValue());
			rm.add(IJ.getImage(), roisG.get(indexFilter.get(j).intValue()),
					sliceList.get(indexFilter.get(j).intValue()).intValue());
		}
		rm.runCommand("Associate", "true");
		rm.runCommand("Centered", "false");
		rm.runCommand("UseNames", "false");
		rm.runCommand(IJ.getImage(), "Show All");
	}

	public void resetFilter() {

		TableRowSorter<TableModel> rowReset = new TableRowSorter<>(modelSlice);
		tableSlice.setRowSorter(rowReset);
		modelS.setValueAt(dataS[0][0], tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(0));
		modelS.setValueAt(dataS[0][1], tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(1));
		modelS.setValueAt(dataS[0][2], tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(2));
		modelS.setValueAt(dataS[0][3], tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(3));
		modelS.setValueAt(dataS[0][4], tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(4));
		modelS.setValueAt(dataS[0][5], tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(5));
		modelS.setValueAt(pValue, tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(6));
		rm.reset();

		for (int i = 0; i < roisG.size(); i++) {
			roisG.get(i).setPosition(sliceList.get(i).intValue());
			rm.add(IJ.getImage(), roisG.get(i), sliceList.get(i).intValue());
		}
		rm.runCommand("Associate", "true");
		rm.runCommand("Centered", "false");
		rm.runCommand("UseNames", "false");
		rm.runCommand(IJ.getImage(), "Show All");

	}

	public void refreshAction() {

		for (int i = 1; i < CellSynapseProcessing_.columnHeadersSlice.length; i++)
			CellSynapseProcessing_.comboSumParam.addItem((String) CellSynapseProcessing_.columnHeadersSlice[i]);

		for (int i = 1; i < columnHeadersSlice.length - 1; i++)
			comboSumParam.addItem((String) columnHeadersSlice[i]);
		prefImages.put(CELLTYPEANALYZER_IMAGES_DEFAULT_PATH, textImages.getText());
		File imageFolder = new File(textImages.getText());
		File[] listOfFiles = imageFolder.listFiles();
		String[] imageTitles = new String[listOfFiles.length];
		imps = new ImagePlus[imageTitles.length];
		icons = new ImageIcon[imps.length];
		for (int i = 0; i < listOfFiles.length; i++) {
			if (listOfFiles[i].isFile())
				imageTitles[i] = listOfFiles[i].getName();

			if (listOfFiles[i].getName().contains(".tif") || listOfFiles[i].getName().contains(".lif")
					|| listOfFiles[i].getName().contains(".lsm") || listOfFiles[i].getName().contains(".czi")) {
				if (listOfFiles[i].getName().contains(".tif") == true
						|| listOfFiles[i].getName().contains(".lsm") == true) {

					imps[i] = IJ.openImage(textImages.getText() + File.separator + imageTitles[i]);
					icons[i] = new ImageIcon(getScaledImage(imps[i].getImage(), 90, 60));
				}

				if (listOfFiles[i].getName().contains(".lif") == true
						|| listOfFiles[i].getName().contains(".czi") == true) {
					lifs = openBF((textImages.getText() + File.separator + imageTitles[i]), false, false, false, false,
							false, true);
					// ImagePlus rgb[] = new ImagePlus[3];

					for (int x = 0; x < lifs.length; x++) {
						// rgb[0] = ChannelSplitter.split(lifs[x])[0];
						// rgb[1] = ChannelSplitter.split(lifs[x])[1];
						// if (ChannelSplitter.split(lifs[x]).length > 2) {
						// rgb[2] = ChannelSplitter.split(lifs[x])[2];
						// IJ.log("passaa");
						// }
						ImagePlus tempLif = null;
						if (ChannelSplitter.split(lifs[x]).length == 2)
							tempLif = new ImagePlus(lifs[x].getTitle(),
									RGBStackMerge.mergeStacks(ChannelSplitter.split(lifs[x])[0].duplicate().getStack(),
											ChannelSplitter.split(lifs[x])[1].duplicate().getStack(), null, false));
						if (ChannelSplitter.split(lifs[x]).length == 3)
							tempLif = new ImagePlus(lifs[x].getTitle(),
									RGBStackMerge.mergeStacks(ChannelSplitter.split(lifs[x])[0].duplicate().getStack(),
											ChannelSplitter.split(lifs[x])[1].duplicate().getStack(),
											ChannelSplitter.split(lifs[x])[2].duplicate().getStack(), false));
						ImageStack[] stacks = null;
						if (ChannelSplitter.split(lifs[x]).length == 4) {
							stacks = new ImageStack[4];
							stacks[0] = ChannelSplitter.split(lifs[x])[0].duplicate().getStack();
							stacks[1] = ChannelSplitter.split(lifs[x])[1].duplicate().getStack();
							stacks[2] = ChannelSplitter.split(lifs[x])[2].duplicate().getStack();
							stacks[3] = ChannelSplitter.split(lifs[x])[3].duplicate().getStack();
							tempLif = new RGBStackMerge().createComposite(lifs[x].getWidth(), lifs[x].getHeight(),
									lifs[x].getNSlices(), stacks, false);
							tempLif.setTitle(lifs[x].getTitle());

						}

						impsLif.add(tempLif);
						iconsLif.add(new ImageIcon(getScaledImage(tempLif.getImage(), 90, 60)));
					}
				}
			}
		}
		if (impsLif.isEmpty() == false) {
			Object[][] dataTImages = new Object[impsLif.size()][columnNames.length];
			for (int i = 0; i < dataTImages.length; i++)
				for (int j = 0; j < dataTImages[i].length; j++)
					dataTImages[i][j] = "";
			modelImages = new DefaultTableModel(dataTImages, columnNames) {

				@Override
				public Class<?> getColumnClass(int column) {
					if (getRowCount() > 0) {
						Object value = getValueAt(0, column);
						if (value != null) {
							return getValueAt(0, column).getClass();
						}
					}

					return super.getColumnClass(column);
				}

				public boolean isCellEditable(int row, int col) {
					return false;
				}

			};
			tableImages.setModel(modelImages);

			for (int i = 0; i < modelImages.getRowCount(); i++) {
				modelImages.setValueAt(iconsLif.get(i), i, tableImages.convertColumnIndexToModel(0));
				modelImages.setValueAt(
						impsLif.get(i).getTitle().substring(0, impsLif.get(i).getTitle().lastIndexOf(".")), i,
						tableImages.convertColumnIndexToModel(1));
				if (impsLif.get(i).getTitle().contains(".tif"))
					modelImages.setValueAt(".tif", i, tableImages.convertColumnIndexToModel(2));
				if (impsLif.get(i).getTitle().contains(".lif"))
					modelImages.setValueAt(".lif", i, tableImages.convertColumnIndexToModel(2));
			}

		}
		if (impsLif.isEmpty() == true) {

			Object[][] dataTImages = new Object[imps.length][columnNames.length];
			for (int i = 0; i < dataTImages.length; i++)
				for (int j = 0; j < dataTImages[i].length; j++)
					dataTImages[i][j] = "";
			modelImages = new DefaultTableModel(dataTImages, columnNames) {

				@Override
				public Class<?> getColumnClass(int column) {
					if (getRowCount() > 0) {
						Object value = getValueAt(0, column);
						if (value != null) {
							return getValueAt(0, column).getClass();
						}
					}

					return super.getColumnClass(column);
				}

				public boolean isCellEditable(int row, int col) {
					return false;
				}

			};
			tableImages.setModel(modelImages);

			for (int i = 0; i < modelImages.getRowCount(); i++) {
				modelImages.setValueAt(icons[i], i, tableImages.convertColumnIndexToModel(0));
				modelImages.setValueAt(imps[i].getTitle().substring(0, imps[i].getTitle().lastIndexOf(".")), i,
						tableImages.convertColumnIndexToModel(1));
				modelImages.setValueAt(imps[i].getTitle().substring(imps[i].getTitle().lastIndexOf(".")), i,
						tableImages.convertColumnIndexToModel(2));
			}

		}
		tableImages.setSelectionBackground(new Color(229, 255, 204));
		tableImages.setSelectionForeground(new Color(0, 102, 0));
		DefaultTableCellRenderer centerRenderer = new DefaultTableCellRenderer();
		centerRenderer.setHorizontalAlignment(JLabel.CENTER);
		tableImages.setDefaultRenderer(String.class, centerRenderer);
		tableImages.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		tableImages.setRowHeight(60);
		tableImages.setAutoCreateRowSorter(true);
		tableImages.getTableHeader().setDefaultRenderer(new SimpleHeaderRenderer());
		tableImages.getColumnModel().getColumn(0).setPreferredWidth(100);
		tableImages.getColumnModel().getColumn(1).setPreferredWidth(450);
		tableImages.getColumnModel().getColumn(2).setPreferredWidth(100);
	}

	public static ImageIcon createImageIcon(String path) {
		java.net.URL imgURL = CellSynapseAnalyzer_.class.getResource(path);
		if (imgURL != null) {
			return new ImageIcon(imgURL);
		} else {
			System.err.println("Couldn't find file: " + path);
			return null;
		}
	}

	public static ImagePlus[] stack2images(ImagePlus imp) {
		String sLabel = imp.getTitle();
		String sImLabel = "";
		ImageStack stack = imp.getStack();

		int sz = stack.getSize();
		int currentSlice = imp.getCurrentSlice(); // to reset ***

		DecimalFormat df = new DecimalFormat("0000"); // for title
		ImagePlus[] arrayOfImages = new ImagePlus[imp.getStack().getSize()];
		for (int n = 1; n <= sz; ++n) {
			imp.setSlice(n); // activate next slice ***

			// Get current image processor from stack. What ever is
			// used here should do a COPY pixels from old processor to
			// new. For instance, ImageProcessor.crop() returns copy.
			ImageProcessor ip = imp.getProcessor(); // ***
			ImageProcessor newip = ip.createProcessor(ip.getWidth(), ip.getHeight());
			newip.setPixels(ip.getPixelsCopy());

			// Create a suitable label, using the slice label if possible
			sImLabel = imp.getStack().getSliceLabel(n);
			if (sImLabel == null || sImLabel.length() < 1) {
				sImLabel = "slice" + df.format(n) + "_" + sLabel;
			}
			// Create new image corresponding to this slice.
			ImagePlus im = new ImagePlus(sImLabel, newip);
			im.setCalibration(imp.getCalibration());
			arrayOfImages[n - 1] = im;

			// Show this image.
			// imp.show();
		}
		// Reset original stack state.
		imp.setSlice(currentSlice); // ***
		if (imp.isProcessor()) {
			ImageProcessor ip = imp.getProcessor();
			ip.setPixels(ip.getPixels()); // ***
		}
		imp.setSlice(currentSlice);
		return arrayOfImages;
	}

	public static double distance(double x1, double y1, double x2, double y2) {
		double x = Math.pow(x2 - x1, 2);
		double y = Math.pow(y2 - y1, 2);
		return Math.sqrt(x + y);
	}

	public static ImagePlus[] openBF(String multiSeriesFileName, boolean splitC, boolean splitT, boolean splitZ,
			boolean autoScale, boolean crop, boolean allSeries) {
		ImporterOptions options;
		ImagePlus[] imps = null;
		try {
			options = new ImporterOptions();
			options.setId(multiSeriesFileName);
			options.setSplitChannels(splitC);
			options.setSplitTimepoints(splitT);
			options.setSplitFocalPlanes(splitZ);
			options.setAutoscale(autoScale);
			options.setStackFormat(ImporterOptions.VIEW_HYPERSTACK);
			options.setStackOrder(ImporterOptions.ORDER_XYCZT);
			options.setCrop(crop);
			options.setOpenAllSeries(allSeries);

			ImportProcess process = new ImportProcess(options);
			if (!process.execute())
				return null;
			DisplayHandler displayHandler = new DisplayHandler(process);
			if (options != null && options.isShowOMEXML()) {
				displayHandler.displayOMEXML();
			}
			List<ImagePlus> impsList = new ImagePlusReaderModified(process).readImages(false);
			imps = impsList.toArray(new ImagePlus[0]);
			if (options != null && options.showROIs()) {
				displayHandler.displayROIs(imps);
			}
			if (!options.isVirtual()) {
				process.getReader().close();
			}

		} catch (Exception e) {

			return null;
		}
		return imps;
	}

}
