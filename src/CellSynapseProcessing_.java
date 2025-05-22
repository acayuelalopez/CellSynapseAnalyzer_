import java.awt.BasicStroke;
import java.awt.Button;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Polygon;
import java.awt.RenderingHints;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.PathIterator;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.prefs.Preferences;
import java.util.regex.PatternSyntaxException;
import java.util.stream.IntStream;

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
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSlider;
import javax.swing.JSpinner;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
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
import org.jfree.chart.ChartPanel;
import org.jfree.chart.plot.IntervalMarker;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.Line;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.gui.ShapeRoi;
import ij.gui.WaitForUserDialog;
import ij.measure.Calibration;
import ij.measure.CurveFitter;
import ij.measure.Measurements;
import ij.plugin.ChannelSplitter;
import ij.plugin.PlugIn;
import ij.plugin.RGBStackMerge;
import ij.plugin.RoiEnlarger;
import ij.plugin.filter.Binary;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import inra.ijpb.morphology.Morphology;
import inra.ijpb.morphology.Reconstruction;
import inra.ijpb.morphology.strel.SquareStrel;
import inra.ijpb.plugins.FillHolesPlugin;
import loci.plugins.in.DisplayHandler;
import loci.plugins.in.ImportProcess;
import loci.plugins.in.ImporterOptions;

public class CellSynapseProcessing_ implements PlugIn {

	JTable tableSlice, tableS;
	DefaultTableModel modelSlice, modelS;
	JScrollPane jScrollPaneSlice;
	ImagePlus[] imps, channels, arrayOfImages, impGSlices, lifs;
	ImageIcon[] icons;
	ImagePlus imp, impCell, impSyn, impProt, impSynMean;
	RoiManager rm;
	Roi[] roisCell, roisCellTotal, roisProt, roisProtTotal, roisSyn, roisSynTotal;
	int indexMax;
	List<Double> sliceList, polList, cellProtList, cellSynList, xCellList, yCellList, xProtList, yProtList, xSynList,
			ySynList, cellArea, proteinArea, synapseArea, aDistList, bDistList, ratioFluorescence;
	List<Double> minDistanceListDouble;
	JCheckBox checkPhysical, checkGraphics;
	JTextField physicalTF;
	JSpinner filterMin, filterMax;
	HistogramFilter hs2;
	IntervalMarker intervalMarker;
	ChartPanel histogram;
	JComboBox comboFilters;
	static JComboBox<String> comboFitters, comboCellCh, comboSynCh, comboProtCh, comboSumParam;
	DefaultListModel<String> modelList;
	JList<String> filterList;
	Double[][] dataSlice, dataS, dataFilter;
	static Object[] columnHeadersS = new String[] { "N", "Sum", "Mean", "Median", "Variance", "stdDev", "Min", "Max",
			"Q1", "Q3", "IQR" },
			columnHeadersSlice = { "Slice-ID", "Pol.Index", "A", "Cell-Protein", "B", "Cell-Synapse", "xCell", "yCell",
					"Cell-Area", "xProtein", "yProtein", "Protein-Area", "xSynapse", "ySynapse", "Synapse-Area",
					"Fluo. Ratio" };
	double pValue;
	double[] xCell, yCell, xProt, yProt, xSyn, ySyn, distUp, distDown;
	List<ImagePlus> impsLif = new ArrayList<ImagePlus>();
	List<ImageIcon> iconsLif = new ArrayList<ImageIcon>();
	JTabbedPane jTabbedPane;
	JPanel panelPicture;
	static JProgressBar pb = new JProgressBar();
	List<Double> distCellProt, distCellSyn, distCellProtN, distCellSynN, aDist, bDist;
	Line[] lineCellProt, lineCellSyn;
	JRadioButton radioCell, radioSyn, radioProt;

	public void showGUI() {

		panelPicture = new JPanel();
		panelPicture.setLayout(new BoxLayout(panelPicture, BoxLayout.Y_AXIS));
		checkPhysical = new JCheckBox("Phy.Units");
		checkGraphics = new JCheckBox("Show Graphics");
		physicalTF = new JTextField(3);
		physicalTF.setEnabled(false);
		JPanel physicalPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		physicalPanel.add(checkPhysical);
		physicalPanel.add(physicalTF);
		physicalPanel.add(checkGraphics);
		JButton buttonProcess = new JButton("");
		ImageIcon iconProcess = createImageIcon("images/processs.png");
		Icon iconProcessCell = new ImageIcon(iconProcess.getImage().getScaledInstance(15, 15, Image.SCALE_SMOOTH));
		buttonProcess.setIcon(iconProcessCell);
		JSeparator separator1 = new JSeparator(SwingConstants.VERTICAL);
		JSeparator separator2 = new JSeparator(SwingConstants.VERTICAL);
		Dimension dime1 = separator1.getPreferredSize();
		Dimension dime2 = separator1.getPreferredSize();
		radioCell = new JRadioButton("Cell-Ch");
		radioCell.setSelected(true);
		radioSyn = new JRadioButton("Synapse-Ch: ");
		radioSyn.setSelected(true);
		radioProt = new JRadioButton("Protein-Ch: ");
		radioProt.setSelected(true);
		comboCellCh = new JComboBox<String>();
		comboCellCh.setPreferredSize(new Dimension(80, 25));
		comboSynCh = new JComboBox<String>();
		comboSynCh.setPreferredSize(new Dimension(80, 25));
		comboProtCh = new JComboBox<String>();
		comboProtCh.setPreferredSize(new Dimension(80, 25));
		JPanel cellPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		cellPanel.add(radioCell);
		cellPanel.add(comboCellCh);
		JPanel synPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		synPanel.add(radioSyn);
		synPanel.add(comboSynCh);
		JPanel protPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		protPanel.add(radioProt);
		protPanel.add(comboProtCh);
		JPanel boxPanel = new JPanel();
		boxPanel.setLayout(new BoxLayout(boxPanel, BoxLayout.Y_AXIS));
		boxPanel.add(cellPanel);
		boxPanel.add(synPanel);
		boxPanel.add(protPanel);
		dime1.height = boxPanel.getPreferredSize().height;
		dime2.height = boxPanel.getPreferredSize().height;
		separator1.setPreferredSize(dime1);
		separator2.setPreferredSize(dime2);
		DefaultTableCellRenderer centerRenderer = new DefaultTableCellRenderer();
		centerRenderer.setHorizontalAlignment(JLabel.CENTER);
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
		JPanel slicePanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		slicePanel.add(jScrollPaneSlice);
		JPanel sPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		sPanel.add(jScrollPaneS);
		JLabel sumParameterLabel = new JLabel("Summary-Parameter: ");
		comboSumParam = new JComboBox<String>();
		JPanel sumParamPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		sumParamPanel.add(sumParameterLabel);
		sumParamPanel.add(comboSumParam);
		JPanel panelBox = new JPanel();
		panelBox.setLayout(new BoxLayout(panelBox, BoxLayout.Y_AXIS));
		panelBox.add(physicalPanel);
		panelBox.add(sumParamPanel);
		// panelPicture.add(imagePanel);
		panelPicture.add(Box.createVerticalStrut(3));
		panelPicture.add(Box.createVerticalStrut(25));
		panelPicture.add(slicePanel);
		panelPicture.add(Box.createVerticalStrut(3));
		physicalPanel.add(buttonProcess);
		JPanel lastPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
		lastPanel.add(separator1);
		lastPanel.add(boxPanel);
		lastPanel.add(separator2);
		lastPanel.add(panelBox);
		panelPicture.add(lastPanel);
		panelPicture.add(sPanel);
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
		ImageIcon iconCS = createImageIcon("images/item1.png");
		Icon CSCell = new ImageIcon(iconCS.getImage().getScaledInstance(16, 16, Image.SCALE_SMOOTH));
		JLabel labelCS = new JLabel("Cell-Synapse Parameters Selection ");
		labelCS.setHorizontalTextPosition(JLabel.TRAILING);
		labelCS.setIcon(CSCell);
		labelCS.setFont(new Font("Verdana", Font.BOLD, 12));
		CellSynapseAnalyzer_.jTabbedPane.addTab(null, CSCell, panelPicture, "Cell-Synapse Analysis");
		CellSynapseAnalyzer_.jTabbedPane.setTabComponentAt(0, labelCS);
		CellSynapseAnalyzer_.jTabbedPane.setMnemonicAt(0, KeyEvent.VK_1);
		CellSynapseAnalyzer_.jTabbedPane.setTabLayoutPolicy(JTabbedPane.SCROLL_TAB_LAYOUT);
		radioCell.addItemListener(new ItemListener() {

			@Override
			public void itemStateChanged(ItemEvent e) {
				if (e.getStateChange() == ItemEvent.DESELECTED)
					comboCellCh.setEnabled(false);
				if (e.getStateChange() == ItemEvent.SELECTED)
					comboCellCh.setEnabled(true);

			}
		});
		radioSyn.addItemListener(new ItemListener() {

			@Override
			public void itemStateChanged(ItemEvent e) {
				if (e.getStateChange() == ItemEvent.DESELECTED)
					comboSynCh.setEnabled(false);
				if (e.getStateChange() == ItemEvent.SELECTED)
					comboSynCh.setEnabled(true);
			}
		});
		radioProt.addItemListener(new ItemListener() {

			@Override
			public void itemStateChanged(ItemEvent e) {
				if (e.getStateChange() == ItemEvent.DESELECTED)
					comboProtCh.setEnabled(false);
				if (e.getStateChange() == ItemEvent.SELECTED)
					comboProtCh.setEnabled(true);
			}
		});
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
				int selectedIndex = comboFilters.getSelectedIndex();

				double values[] = new double[dataSlice.length];
				for (int r = 0; r < dataSlice.length; r++)
					values[r] = dataSlice[r][selectedIndex];

				double maxValue = values[0];
				double minValue = values[0];
				for (int i = 1; i < values.length; i++) {
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

		buttonProcess.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				processTool();

			}
		});
		btnCF.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				curveFitter();
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
				FileWriter csv = new FileWriter(new File(fileToSave.getAbsolutePath() + File.separator
						+ CellSynapseAnalyzer_.imageLabel.getText() + "_data.csv"));

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

		List<Double> listFiltered = new ArrayList<Double>();
		for (int row = 0; row < tableSlice.getRowCount(); row++)
			listFiltered.add((Double) tableSlice.getModel().getValueAt(tableSlice.convertRowIndexToModel(row),
					comboSumParam.getSelectedIndex() + 1));

		modelS.setValueAt((double) tableSlice.getRowSorter().getViewRowCount(), tableS.convertRowIndexToModel(0),
				tableS.convertColumnIndexToModel(0));
		modelS.setValueAt(Math.round(listFiltered.stream().mapToDouble(Double::doubleValue).sum() * 1000.0) / 1000.0,
				tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(1));
		modelS.setValueAt(
				Math.round((listFiltered.stream().mapToDouble(Double::doubleValue).sum()) / polList.size() * 1000.0)
						/ 1000.0,
				tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(2));
		modelS.setValueAt(Math.round(median(listFiltered) * 1000.0) / 1000.0, tableS.convertRowIndexToModel(0),
				tableS.convertColumnIndexToModel(3));
		modelS.setValueAt(Math.round(variance(listFiltered) * 1000.0) / 1000.0, tableS.convertRowIndexToModel(0),
				tableS.convertColumnIndexToModel(4));
		modelS.setValueAt(Math.round(sd(listFiltered) * 1000.0) / 1000.0, tableS.convertRowIndexToModel(0),
				tableS.convertColumnIndexToModel(5));
		modelS.setValueAt(Math.round(Collections.min(listFiltered) * 1000.0) / 1000.0, tableS.convertRowIndexToModel(0),
				tableS.convertColumnIndexToModel(6));
		modelS.setValueAt(Math.round(Collections.max(listFiltered) * 1000.0) / 1000.0, tableS.convertRowIndexToModel(0),
				tableS.convertColumnIndexToModel(7));
		Collections.sort(listFiltered);
		int mid_index = medianIQR(listFiltered, 0, listFiltered.size());
		double Q1 = listFiltered.get(medianIQR(listFiltered, 0, mid_index));
		double Q3 = listFiltered.get(medianIQR(listFiltered, mid_index + 1, listFiltered.size()));
		double IQR = (Q3 - Q1);
		modelS.setValueAt(Math.round(Q1 * 1000.0) / 1000.0, tableS.convertRowIndexToModel(0),
				tableS.convertColumnIndexToModel(8));
		modelS.setValueAt(Math.round(Q3 * 1000.0) / 1000.0, tableS.convertRowIndexToModel(0),
				tableS.convertColumnIndexToModel(9));
		modelS.setValueAt(Math.round(IQR * 1000.0) / 1000.0, tableS.convertRowIndexToModel(0),
				tableS.convertColumnIndexToModel(10));
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
		modelS.setValueAt(dataS[0][6], tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(6));
		modelS.setValueAt(dataS[0][7], tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(7));
		modelS.setValueAt(dataS[0][8], tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(8));
		modelS.setValueAt(dataS[0][9], tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(9));
		modelS.setValueAt(dataS[0][10], tableS.convertRowIndexToModel(0), tableS.convertColumnIndexToModel(10));

	}

	public void processTool() {
		Thread process0 = new Thread(new Runnable() {

			public void run() {
				rm = RoiManager.getInstance();
				if (null == rm)
					rm = new RoiManager();
				if (null != rm)
					rm.close();
				List<ImagePlus> listImages = new ArrayList<ImagePlus>();
				sliceList = new ArrayList<Double>();
				polList = new ArrayList<Double>();
				aDistList = new ArrayList<Double>();
				bDistList = new ArrayList<Double>();
				cellProtList = new ArrayList<Double>();
				cellSynList = new ArrayList<Double>();
				xCellList = new ArrayList<Double>();
				yCellList = new ArrayList<Double>();
				xProtList = new ArrayList<Double>();
				yProtList = new ArrayList<Double>();
				xSynList = new ArrayList<Double>();
				ySynList = new ArrayList<Double>();
				cellArea = new ArrayList<Double>();
				proteinArea = new ArrayList<Double>();
				synapseArea = new ArrayList<Double>();
				ratioFluorescence = new ArrayList<Double>();
				minDistanceListDouble = new ArrayList<Double>();

				imp = CellSynapseAnalyzer_.imp;
				if (imp == null) {
					IJ.error("You should select and  an image from your directory.");
					return;
				}
				if (imp != null) {
					channels = ChannelSplitter.split(imp);
					impCell = channels[comboCellCh.getSelectedIndex()].duplicate();
					impSyn = channels[comboSynCh.getSelectedIndex()].duplicate();
					impSynMean = channels[comboSynCh.getSelectedIndex()].duplicate();
					impProt = channels[comboProtCh.getSelectedIndex()].duplicate();
					listImages.add(imp.duplicate());
					for (int i = 0; i < channels.length; i++)
						listImages.add(channels[i]);
					arrayOfImages = new ImagePlus[listImages.size()];
					for (int i = 0; i < listImages.size(); i++)
						arrayOfImages[i] = listImages.get(i);
					listImages.toArray(arrayOfImages);
					IJ.setTool("rectangle");
					WaitForUserDialog cellDialog = new WaitForUserDialog("Draw cell area, then click OK.");
					WaitForUserDialog protDialog = new WaitForUserDialog("Draw protein area, then click OK.");
					WaitForUserDialog synDialog = new WaitForUserDialog("Draw synapse area, then click OK.");

					if (radioCell.isSelected() == Boolean.TRUE) {
						cellDialog.show();
						cellAction(impCell);
					}
					IJ.setTool("rectangle");
					if (radioProt.isSelected() == Boolean.TRUE) {
						protDialog.show();
						proteinAction(impProt);
					}
					if (radioSyn.isSelected() == Boolean.TRUE) {
						synDialog.show();
						synapseAction(impSyn);
					}
					processingAction();
					pb.setString("Done!!");

				}
			}
		});
		process0.start();
	}

	public void processingAction() {

		// Cell Processing
		distCellProt = new ArrayList<Double>();
		distCellSyn = new ArrayList<Double>();
		distCellProtN = new ArrayList<Double>();
		distCellSynN = new ArrayList<Double>();
		aDist = new ArrayList<Double>();
		bDist = new ArrayList<Double>();
		// IJ.run(impCell, "Auto Threshold", "method=Huang2 ignore_black white stack");
		// IJ.run(impCell, "Close-", "stack");
		ImagePlus[] impCellSlices = CellSynapseAnalyzer_.stack2images(impCell);
		// roisCell = new Roi[impCellSlices.length];
		// roisCellTotal = new Roi[impCellSlices.length];
		// int[] indexesMax = new int[roisCellTotal.length];
		xCell = new double[roisCellTotal.length];
		yCell = new double[roisCellTotal.length];
		xProt = new double[roisProtTotal.length];
		yProt = new double[roisProtTotal.length];
		xSyn = new double[roisSynTotal.length];
		ySyn = new double[roisSynTotal.length];
		lineCellProt = new Line[roisCellTotal.length];
		lineCellSyn = new Line[roisCellTotal.length];
		distUp = new double[roisCellTotal.length];
		distDown = new double[roisCellTotal.length];
		final int MAX = impCellSlices.length;
		pb.setMinimum(0);
		pb.setMaximum(MAX);
		pb.setStringPainted(true);
		for (int i = 0; i < impCellSlices.length; i++) {
			List<Double> xSynTotal = new ArrayList<Double>();
			List<Double> ySynTotal = new ArrayList<Double>();
			List<Double> x3 = new ArrayList<Double>();
			List<Double> y3 = new ArrayList<Double>();
			List<Double> xDown = new ArrayList<Double>();
			List<Double> yDown = new ArrayList<Double>();
			List<Double> xUp = new ArrayList<Double>();
			List<Double> yUp = new ArrayList<Double>();

			final int currentValue = i + 1;
			try {
				SwingUtilities.invokeLater(new Runnable() {
					public void run() {
						pb.setString("Processing... " + (currentValue * 100) / impCellSlices.length + "%");
						pb.setValue(currentValue);

					}
				});
				java.lang.Thread.sleep(100);
			} catch (InterruptedException e1) {
			}
			xCell[i] = roisCellTotal[i].getStatistics().xCentroid;
			yCell[i] = roisCellTotal[i].getStatistics().yCentroid;
			xProt[i] = roisProtTotal[i].getStatistics().xCentroid;
			yProt[i] = roisProtTotal[i].getStatistics().yCentroid;

			// Cell-Protein Distance
			lineCellProt[i] = new Line(xCell[i], yCell[i], xProt[i], yProt[i]);
			if (yProt[i] < yCell[i]) {
				// IJ.log("pasa menos 1");
				distCellProt.add(Math.sqrt(((xCell[i] - xProt[i]) * (xCell[i] - xProt[i]))
						+ ((yCell[i] - yProt[i]) * (yCell[i] - yProt[i]))));

				x3.add(x3Coordinate(xProt[i], xCell[i], yProt[i], yCell[i], 1000.0));
				y3.add(y3Coordinate(xProt[i], xCell[i], yProt[i], yCell[i], 1000.0));
				Roi roiIntersection = (new ShapeRoi(roisSynTotal[i]).and(new ShapeRoi(
						new Line(xCell[i], yCell[i], x3Coordinate(xProt[i], xCell[i], yProt[i], yCell[i], 1000.0),
								y3Coordinate(xProt[i], xCell[i], yProt[i], yCell[i], 1000.0))))).shapeToRoi();
				if (roiIntersection == null) {
					roiIntersection = (new ShapeRoi(roisSyn[i]).and(new ShapeRoi(
							new Line(xCell[i], yCell[i], x3Coordinate(xProt[i], xCell[i], yProt[i], yCell[i], 1000.0),
									y3Coordinate(xProt[i], xCell[i], yProt[i], yCell[i], 1000.0))))).shapeToRoi();
					roisSynTotal[i] = roisSyn[i];
				}

				for (int x = 0; x < roiIntersection.getFloatPolygon().xpoints.length; x++)
					if (roisSynTotal[i].containsPoint(roiIntersection.getFloatPolygon().xpoints[x],
							roiIntersection.getFloatPolygon().ypoints[x]) == true) {
						xSynTotal.add((double) roiIntersection.getFloatPolygon().xpoints[x]);
						ySynTotal.add((double) roiIntersection.getFloatPolygon().ypoints[x]);

					}

				xSyn[i] = xSynTotal.get(0);
				ySyn[i] = ySynTotal.get(0);
				lineCellSyn[i] = new Line(xCell[i], yCell[i], xSyn[i], ySyn[i]);
				distCellSyn.add(Math.sqrt(
						((xSyn[i] - xCell[i]) * (xSyn[i] - xCell[i])) + ((ySyn[i] - yCell[i]) * (ySyn[i] - yCell[i]))));
				// Up
				Roi roiIntersectionUp = (new ShapeRoi(roisCellTotal[i]).and(new ShapeRoi(
						new Line(xCell[i], yCell[i], x3Coordinate(xCell[i], xProt[i], yCell[i], yProt[i], 1000.0),
								y3Coordinate(xCell[i], xProt[i], yCell[i], yProt[i], 1000.0))))).shapeToRoi();
				for (int x = 0; x < roiIntersectionUp.getFloatPolygon().xpoints.length; x++) {
					xUp.add((double) roiIntersectionUp.getFloatPolygon().xpoints[x]);
					yUp.add((double) roiIntersectionUp.getFloatPolygon().ypoints[x]);
				}

				int minIndex = yUp.indexOf(Collections.min(yUp));
				distUp[i] = Math.sqrt(((xUp.get(minIndex) - xCell[i]) * (xUp.get(minIndex) - xCell[i]))
						+ ((yUp.get(minIndex) - yCell[i]) * (yUp.get(minIndex) - yCell[i])));
				distCellProtN.add((-1.0) * (distCellProt.get(i) / distUp[i]));
				// IJ.log(distUp[i] + "-----distUp---" + i + "---" + distCellProt.get(i));
				// Down
				Roi roiIntersectionDown = (new ShapeRoi(roisCellTotal[i]).and(new ShapeRoi(
						new Line(xCell[i], yCell[i], x3Coordinate(xCell[i], xSyn[i], yCell[i], ySyn[i], 1000.0),
								y3Coordinate(xCell[i], xSyn[i], yCell[i], ySyn[i], 1000.0))))).shapeToRoi();
				for (int x = 0; x < roiIntersectionDown.getFloatPolygon().xpoints.length; x++) {
					xDown.add((double) roiIntersectionDown.getFloatPolygon().xpoints[x]);
					yDown.add((double) roiIntersectionDown.getFloatPolygon().ypoints[x]);
				}

				int maxIndex = yDown.indexOf(Collections.max(yDown));
				distDown[i] = Math.sqrt(((xDown.get(maxIndex) - xCell[i]) * (xDown.get(maxIndex) - xCell[i]))
						+ ((yDown.get(maxIndex) - yCell[i]) * (yDown.get(maxIndex) - yCell[i])));
				distCellSynN.add(distCellSyn.get(i) / distDown[i]);
				aDist.add(((-1.0) * (distCellProt.get(i) / distUp[i])));
				bDist.add(distCellSyn.get(i) / distDown[i]);

			}
			if (yProt[i] >= yCell[i]) {
				// IJ.log("pasa normal");
				distCellProt.add(Math.sqrt(((xProt[i] - xCell[i]) * (xProt[i] - xCell[i]))
						+ ((yProt[i] - yCell[i]) * (yProt[i] - yCell[i]))));

				x3.add(x3Coordinate(xCell[i], xProt[i], yCell[i], yProt[i], 1000.0));
				y3.add(y3Coordinate(xCell[i], xProt[i], yCell[i], yProt[i], 1000.0));
				Roi roiIntersection = (new ShapeRoi(roisSynTotal[i]).and(new ShapeRoi(
						new Line(xCell[i], yCell[i], x3Coordinate(xCell[i], xProt[i], yCell[i], yProt[i], 100000.0),
								y3Coordinate(xCell[i], xProt[i], yCell[i], yProt[i], 100000.0))))).shapeToRoi();
				if (roiIntersection == null) {
					roiIntersection = (new ShapeRoi(roisSyn[i]).and(new ShapeRoi(
							new Line(xCell[i], yCell[i], x3Coordinate(xCell[i], xProt[i], yCell[i], yProt[i], 100000.0),
									y3Coordinate(xCell[i], xProt[i], yCell[i], yProt[i], 100000.0))))).shapeToRoi();

				}
				if (roiIntersection != null) {

					for (int x = 0; x < roiIntersection.getFloatPolygon().xpoints.length; x++)
						ySynTotal.add((double) roiIntersection.getFloatPolygon().ypoints[x]);

					xSyn[i] = roiIntersection.getFloatPolygon().xpoints[ySynTotal.indexOf(Collections.min(ySynTotal))];
					ySyn[i] = roiIntersection.getFloatPolygon().ypoints[ySynTotal.indexOf(Collections.min(ySynTotal))];

					lineCellSyn[i] = new Line(xCell[i], yCell[i], xSyn[i], ySyn[i]);
					distCellSyn.add(Math.sqrt(((xSyn[i] - xCell[i]) * (xSyn[i] - xCell[i]))
							+ ((ySyn[i] - yCell[i]) * (ySyn[i] - yCell[i]))));
					// Down
					Roi roiIntersectionDown = (new ShapeRoi(roisCellTotal[i]).and(new ShapeRoi(
							new Line(xCell[i], yCell[i], x3Coordinate(xCell[i], xProt[i], yCell[i], yProt[i], 1000.0),
									y3Coordinate(xCell[i], xProt[i], yCell[i], yProt[i], 1000.0))))).shapeToRoi();
					for (int x = 0; x < roiIntersectionDown.getFloatPolygon().xpoints.length; x++) {
						xDown.add((double) roiIntersectionDown.getFloatPolygon().xpoints[x]);
						yDown.add((double) roiIntersectionDown.getFloatPolygon().ypoints[x]);
					}

					int maxIndex = yDown.indexOf(Collections.max(yDown));
					distDown[i] = Math.sqrt(((xDown.get(maxIndex) - xCell[i]) * (xDown.get(maxIndex) - xCell[i]))
							+ ((yDown.get(maxIndex) - yCell[i]) * (yDown.get(maxIndex) - yCell[i])));
					distCellProtN.add(distCellProt.get(i) / distDown[i]);
					distCellSynN.add(distCellSyn.get(i) / distDown[i]);
					aDist.add(distCellProt.get(i) / distDown[i]);
					bDist.add(distCellSyn.get(i) / distDown[i]);
				}

			}

			sliceList.add((double) i + 1.0);
			if (checkPhysical.isSelected() == true) {
				double pixelSize = Double.parseDouble(physicalTF.getText());
				polList.add(Math.round(pixelSize * (aDist.get(i) / bDist.get(i)) * 1000.0) / 1000.0);
				aDistList.add(Math.round(pixelSize * (aDist.get(i)) * 1000.0) / 1000.0);
				bDistList.add(Math.round(pixelSize * (bDist.get(i)) * 1000.0) / 1000.0);
				cellProtList.add(Math.round(pixelSize * distCellProt.get(i) * 1000.0) / 1000.0);
				cellSynList.add(Math.round(pixelSize * distCellSyn.get(i) * 1000.0) / 1000.0);
				xCellList.add(Math.round(pixelSize * xCell[i] * 1000.0) / 1000.0);
				yCellList.add(Math.round(pixelSize * yCell[i] * 1000.0) / 1000.0);
				xProtList.add(Math.round(pixelSize * xProt[i] * 1000.0) / 1000.0);
				yProtList.add(Math.round(pixelSize * yProt[i] * 1000.0) / 1000.0);
				xSynList.add(Math.round(pixelSize * xSyn[i] * 1000.0) / 1000.0);
				ySynList.add(Math.round(pixelSize * ySyn[i] * 1000.0) / 1000.0);
				cellArea.add(Math.round(pixelSize * roisCellTotal[i].getStatistics().area * 1000.0) / 1000.0);
				proteinArea.add(Math.round(pixelSize * roisProtTotal[i].getStatistics().area * 1000.0) / 1000.0);
				synapseArea.add(Math.round(pixelSize * roisSynTotal[i].getStatistics().area * 1000.0) / 1000.0);

			}
			ImagePlus[] impSynMeanSlices = CellSynapseAnalyzer_.stack2images(impSynMean);
			if (checkPhysical.isSelected() == false) {
				// IJ.log((distCellProtN.get(i) + "---" + distCellSynN.get(i)));

				polList.add(Math.round((aDist.get(i) / bDist.get(i)) * 1000.0) / 1000.0);
				aDistList.add(Math.round((aDist.get(i)) * 1000.0) / 1000.0);
				bDistList.add(Math.round((bDist.get(i)) * 1000.0) / 1000.0);
				cellProtList.add(Math.round(distCellProt.get(i) * 1000.0) / 1000.0);
				cellSynList.add(Math.round(distCellSyn.get(i) * 1000.0) / 1000.0);
				xCellList.add(Math.round(xCell[i] * 1000.0) / 1000.0);
				yCellList.add(Math.round(yCell[i] * 1000.0) / 1000.0);
				xProtList.add(Math.round(xProt[i] * 1000.0) / 1000.0);
				yProtList.add(Math.round(yProt[i] * 1000.0) / 1000.0);
				xSynList.add(Math.round(xSyn[i] * 1000.0) / 1000.0);
				ySynList.add(Math.round(ySyn[i] * 1000.0) / 1000.0);
				cellArea.add(Math.round(roisCellTotal[i].getStatistics().area * 1000.0) / 1000.0);
				proteinArea.add(Math.round(roisProtTotal[i].getStatistics().area * 1000.0) / 1000.0);
				synapseArea.add(Math.round(roisSynTotal[i].getStatistics().area * 1000.0) / 1000.0);
			}
			IJ.run(impSynMeanSlices[i], "8-bit", "");
			impSynMeanSlices[i].setRoi(roisSynTotal[i]);
			impSynMeanSlices[i].setRoi(roisCellTotal[i]);
			ratioFluorescence.add(
					Math.round((roisSynTotal[i].getStatistics().mean / roisCellTotal[i].getStatistics().mean) * 1000.0)
							/ 1000.0);
			// IJ.log(roisSynTotal[i].getStatistics().mean + "-----" +
			// roisCellTotal[i].getStatistics().mean);

		}
		ArrayList<List<Double>> listOfLists = new ArrayList<List<Double>>();
		dataSlice = new Double[xCellList.size()][16];
		for (int i = 0; i < dataSlice.length; i++) {
			dataSlice[i][0] = sliceList.get(i);
			dataSlice[i][1] = polList.get(i);
			dataSlice[i][2] = aDistList.get(i);
			dataSlice[i][3] = cellProtList.get(i);
			dataSlice[i][4] = bDistList.get(i);
			dataSlice[i][5] = cellSynList.get(i);
			dataSlice[i][6] = xCellList.get(i);
			dataSlice[i][7] = yCellList.get(i);
			dataSlice[i][8] = cellArea.get(i);
			dataSlice[i][9] = xProtList.get(i);
			dataSlice[i][10] = yProtList.get(i);
			dataSlice[i][11] = proteinArea.get(i);
			dataSlice[i][12] = xSynList.get(i);
			dataSlice[i][13] = ySynList.get(i);
			dataSlice[i][14] = synapseArea.get(i);
			dataSlice[i][15] = ratioFluorescence.get(i);

		}
		listOfLists.add(polList);
		listOfLists.add(aDistList);
		listOfLists.add(cellProtList);
		listOfLists.add(bDistList);
		listOfLists.add(cellSynList);
		listOfLists.add(xCellList);
		listOfLists.add(yCellList);
		listOfLists.add(cellArea);
		listOfLists.add(xProtList);
		listOfLists.add(yProtList);
		listOfLists.add(proteinArea);
		listOfLists.add(xSynList);
		listOfLists.add(ySynList);
		listOfLists.add(synapseArea);
		listOfLists.add(ratioFluorescence);

		modelSlice = new DefaultTableModel(dataSlice, columnHeadersSlice) {

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

		DefaultTableCellRenderer centerRenderer = new DefaultTableCellRenderer();
		centerRenderer.setHorizontalAlignment(JLabel.CENTER);
		tableSlice.setDefaultRenderer(Double.class, centerRenderer);
		tableSlice.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		tableSlice.setRowHeight(40);
		tableSlice.getTableHeader().setDefaultRenderer(new SimpleHeaderRenderer());
		for (int u = 0; u < tableSlice.getColumnCount(); u++)
			tableSlice.getColumnModel().getColumn(u).setPreferredWidth(130);
		if (checkGraphics.isSelected() == Boolean.TRUE) {
			Overlay overlay = new Overlay();
			for (int s = 0; s < impCellSlices.length; s++) {

				roisCellTotal[s].setPosition(0, s + 1, 0);
				roisCellTotal[s].setStrokeColor(Color.RED);
				overlay.add(roisCellTotal[s], "Cell");
				roisProtTotal[s].setPosition(0, s + 1, 0);
				roisProtTotal[s].setStrokeColor(Color.PINK);
				overlay.add(roisProtTotal[s], "Protein");
				roisSynTotal[s].setPosition(0, s + 1, 0);
				roisSynTotal[s].setStrokeColor(Color.BLUE);
				overlay.add(roisSynTotal[s], "Synapse");
				lineCellProt[s].setPosition(0, s + 1, 0);
				lineCellProt[s].setStrokeColor(Color.LIGHT_GRAY);
				overlay.add(lineCellProt[s], "Cell-Prot: A");
				lineCellSyn[s].setPosition(0, s + 1, 0);
				lineCellSyn[s].setStrokeColor(Color.ORANGE);
				overlay.add(lineCellSyn[s], "Cell-Syn: B");
			}
			imp.setOverlay(overlay);
			imp.killRoi();
		}

		dataS = new Double[1][11];
		// Collections.sort(minDistanceListDouble);
		dataS[0][0] = (double) sliceList.size();

		dataS[0][1] = Math
				.round(listOfLists.get(comboSumParam.getSelectedIndex()).stream().mapToDouble(Double::doubleValue).sum()
						* 1000.0)
				/ 1000.0;
		dataS[0][2] = Math.round(
				(listOfLists.get(comboSumParam.getSelectedIndex()).stream().mapToDouble(Double::doubleValue).sum())
						/ listOfLists.get(comboSumParam.getSelectedIndex()).size() * 1000.0)
				/ 1000.0;
		dataS[0][3] = Math.round(median(listOfLists.get(comboSumParam.getSelectedIndex())) * 1000.0) / 1000.0;
		dataS[0][4] = Math.round(variance(listOfLists.get(comboSumParam.getSelectedIndex())) * 1000.0) / 1000.0;
		dataS[0][5] = Math.round(sd(listOfLists.get(comboSumParam.getSelectedIndex())) * 1000.0) / 1000.0;
		dataS[0][6] = Math.round(Collections.min(listOfLists.get(comboSumParam.getSelectedIndex())) * 1000.0) / 1000.0;
		dataS[0][7] = Math.round(Collections.max(listOfLists.get(comboSumParam.getSelectedIndex())) * 1000.0) / 1000.0;
		dataS[0][8] = Math.round(Q1(listOfLists.get(comboSumParam.getSelectedIndex()),
				listOfLists.get(comboSumParam.getSelectedIndex()).size()) * 1000.0) / 1000.0;
		dataS[0][9] = Math.round(Q3(listOfLists.get(comboSumParam.getSelectedIndex()),
				listOfLists.get(comboSumParam.getSelectedIndex()).size()) * 1000.0) / 1000.0;
		dataS[0][10] = Math.round(IQR(listOfLists.get(comboSumParam.getSelectedIndex()),
				listOfLists.get(comboSumParam.getSelectedIndex()).size()) * 1000.0) / 1000.0;

		modelS = new DefaultTableModel(dataS, columnHeadersS) {

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
		tableS.setModel(modelS);
		tableS.setSelectionBackground(new Color(229, 255, 204));
		tableS.setSelectionForeground(new Color(0, 102, 0));
		tableS.setDefaultRenderer(Double.class, centerRenderer);
		tableS.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		tableS.setRowHeight(40);
		tableS.getTableHeader().setDefaultRenderer(new SimpleHeaderRenderer());
		for (int u = 0; u < tableS.getColumnCount(); u++)
			tableS.getColumnModel().getColumn(u).setPreferredWidth(130);
		IJ.run("Labels...", "color=white font=12 show use");
	}

	public void synapseAction(ImagePlus impSyn) {

		Roi synRoi = imp.getRoi();
		// Synapse Protein
		IJ.run(impSyn, "Auto Threshold", "method=MaxEntropy ignore_black white stack");
		ImagePlus[] impSynSlices = CellSynapseAnalyzer_.stack2images(impSyn);
		roisSyn = new Roi[impSynSlices.length];
		Roi[] roisSynInitial = new Roi[impSynSlices.length];
		roisSynTotal = new Roi[impSynSlices.length];
		int[] indexesMax = new int[impSynSlices.length];
		for (int i = 0; i < impSynSlices.length; i++) {
			impSynSlices[i] = new ImagePlus(impSynSlices[i].getTitle(),
					Morphology.closing(impSynSlices[i].getProcessor(), SquareStrel.fromDiameter(2)));
			IJ.run(impSynSlices[i], "Create Selection", "");
			roisSyn[i] = impSynSlices[i].getRoi();
			roisSynInitial[i] = RoiEnlarger.enlarge(
					new ShapeRoi(roisSyn[i]).xor((new ShapeRoi(synRoi).or(new ShapeRoi(roisSyn[i])))).shapeToRoi(),
					2.0);

			Roi[] roisSynSplit = new ShapeRoi(roisSynInitial[i]).getRois();

			double[] roisAreaSyn = new double[roisSynSplit.length];
			for (int j = 0; j < roisSynSplit.length; j++)
				roisAreaSyn[j] = roisSynSplit[j].getStatistics().area;

			double max = roisAreaSyn[0];
			indexMax = 0;

			for (int j = 0; j < roisAreaSyn.length; j++) {
				if (max < roisAreaSyn[j]) {
					max = roisAreaSyn[j];
					indexMax = j;
					indexesMax[i] = indexMax;

				}
				//if (roisAreaSyn[j] >= roisAreaSyn[indexesMax[i]] / 3)
					//roisSynTotal[i] = new ShapeRoi(roisSynSplit[j]).or(new ShapeRoi(roisSynSplit[indexesMax[i]]))
							//.shapeToRoi();
				roisSynTotal[i] = roisSynInitial[i];
			}

		}
	//	for (int i = 0; i < roisCellTotal.length; i++)
			//roisCellTotal[i] = new ShapeRoi(roisCellTotal[i]).xor((new ShapeRoi(roisCellTotal[i]).and(new ShapeRoi(roisSynTotal[i]))))
					//.shapeToRoi();

	}

	public void proteinAction(ImagePlus impProt) {
		Roi protRoi = imp.getRoi();
		// Protein Processing
		IJ.run(impProt, "Auto Threshold", "method=Otsu ignore_black white stack");
		ImagePlus[] impProtSlices = CellSynapseAnalyzer_.stack2images(impProt);
		roisProt = new Roi[impProtSlices.length];
		roisProtTotal = new Roi[impProtSlices.length];
		for (int i = 0; i < impProtSlices.length; i++) {
			IJ.run(impProtSlices[i], "Create Selection", "");
			roisProt[i] = impProtSlices[i].getRoi();
			roisProtTotal[i] = new ShapeRoi(roisProt[i]).xor((new ShapeRoi(protRoi).or(new ShapeRoi(roisProt[i]))))
					.shapeToRoi();

		}

	}

	public void cellAction(ImagePlus impCell) {
		Roi cellRoi = imp.getRoi();
		// Cell Processing
		IJ.run(impCell, "Auto Threshold", "method=Huang ignore_black white stack");
		//IJ.run(impCell, "Erode", "stack");
		ImagePlus[] impCellSlices = CellSynapseAnalyzer_.stack2images(impCell);
		roisCell = new Roi[impCellSlices.length];
		roisCellTotal = new Roi[impCellSlices.length];
		int[] indexesMax = new int[roisCellTotal.length];
		for (int i = 0; i < impCellSlices.length; i++) {
			IJ.run(impCellSlices[i], "Create Selection", "");
			roisCell[i] = impCellSlices[i].getRoi();
			roisCellTotal[i] = new ShapeRoi(cellRoi).xor((new ShapeRoi(cellRoi).and(new ShapeRoi(roisCell[i]))))
					.shapeToRoi();
			Roi[] roisCellSplit = new ShapeRoi(roisCellTotal[i]).getRois();
			List<Roi> roisCellSplitList = new ArrayList<Roi>();

			for (int j = 0; j < roisCellSplit.length; j++)
				if (roisCellSplit[j].getType() != Roi.RECTANGLE)
					roisCellSplitList.add(roisCellSplit[j]);
			double[] roisAreaCell = new double[roisCellSplitList.size()];
			for (int x = 0; x < roisCellSplitList.size(); x++)
				roisAreaCell[x] = roisCellSplitList.get(x).getStatistics().area;

			double max = roisAreaCell[0];
			indexMax = 0;
			for (int j = 0; j < roisAreaCell.length; j++) {
				if (max < roisAreaCell[j]) {
					max = roisAreaCell[j];
					indexMax = j;
					indexesMax[i] = indexMax;

				}
			}

			roisCellTotal[i] = roisCellSplitList.get(indexesMax[i]);

		}

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

	public static Double findMin(List<Double> list) {

		// check list is empty or not
		if (list == null || list.size() == 0) {
			return Double.MAX_VALUE;
		}

		List<Double> sortedlist = new ArrayList<>(list);
		Collections.sort(sortedlist);
		return sortedlist.get(0);
	}

	public static double distance(double x1, double y1, double x2, double y2) {
		double x = Math.pow(x2 - x1, 2);
		double y = Math.pow(y2 - y1, 2);
		return Math.sqrt(x + y);
	}

	public static double median(List<Double> a) {
		int middle = a.size() / 2;

		if (a.size() % 2 == 1) {
			return a.get(middle);
		} else {
			return (a.get(middle - 1) + a.get(middle)) / 2.0;
		}
	}

	public static double IQR(List<Double> a, int n) {

		Collections.sort(a);
		int mid_index = medianIQR(a, 0, n);
		double Q1 = a.get(medianIQR(a, 0, mid_index));
		double Q3 = a.get(medianIQR(a, mid_index + 1, n));

		return (Q3 - Q1);
	}

	public static double Q3(List<Double> a, int n) {
		Collections.sort(a);
		int mid_index = medianIQR(a, 0, n);
		double Q3 = a.get(medianIQR(a, mid_index + 1, n));
		return (Q3);
	}

	public static double Q1(List<Double> a, int n) {
		Collections.sort(a);
		int mid_index = medianIQR(a, 0, n);
		double Q1 = a.get(medianIQR(a, 0, mid_index));
		return (Q1);
	}

	public static double variance(List<Double> a) {
		double mean = a.stream().mapToDouble(Double::doubleValue).sum() / a.size();
		double square = 0;
		for (double element : a)
			square += (element - mean) * (element - mean);
		return square / a.size();
	}

	public static int medianIQR(List<Double> a, int d, int mid_index) {
		int n = mid_index - d + 1;
		n = (n + 1) / 2 - 1;
		return n + d;
	}

	public static double sd(List<Double> a) {
		int sum = 0;
		double mean = a.stream().mapToDouble(Double::doubleValue).sum() / a.size();

		for (Double i : a)
			sum += Math.pow((i - mean), 2);
		return Math.sqrt(sum / (a.size() - 1)); // sample
	}

	public int smallestIndex(float[] array) {
		float currentValue = array[0];
		int smallestIndex = 0;
		for (int j = 1; j < array.length; j++) {
			if (array[j] < currentValue) {
				currentValue = array[j];
				smallestIndex = j;
			}
		}

		return smallestIndex;
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

	public double x3Coordinate(double xCell, double xProt, double yCell, double yProt, double n) {
		double d = Math.sqrt(((xProt - xCell) * (xProt - xCell)) + ((yProt - yCell) * (yProt - yCell))); // distance
		double r = n / d; // segment ratio
		double x3 = r * xProt + (1 - r) * xCell; // find point that divides the segment
		return x3;
	}

	public double y3Coordinate(double xCell, double xProt, double yCell, double yProt, double n) {
		double d = Math.sqrt(((xProt - xCell) * (xProt - xCell)) + ((yProt - yCell) * (yProt - yCell))); // distance
		double r = n / d; // segment ratio
		double y3 = r * yProt + (1 - r) * yCell; // into the ratio (1-r):r
		return y3;
	}

	public static Set<Point2D> getIntersections(final Polygon poly, final Line2D.Double line) throws Exception {

		final PathIterator polyIt = poly.getPathIterator(null); // Getting an iterator along the polygon path
		final double[] coords = new double[6]; // Double array with length 6 needed by iterator
		final double[] firstCoords = new double[2]; // First point (needed for closing polygon path)
		final double[] lastCoords = new double[2]; // Previously visited point
		final Set<Point2D> intersections = new HashSet<Point2D>(); // List to hold found intersections
		polyIt.currentSegment(firstCoords); // Getting the first coordinate pair
		lastCoords[0] = firstCoords[0]; // Priming the previous coordinate pair
		lastCoords[1] = firstCoords[1];
		polyIt.next();
		while (!polyIt.isDone()) {
			final int type = polyIt.currentSegment(coords);
			switch (type) {
			case PathIterator.SEG_LINETO: {
				final Line2D.Double currentLine = new Line2D.Double(lastCoords[0], lastCoords[1], coords[0], coords[1]);
				if (currentLine.intersectsLine(line))
					intersections.add(getIntersection(currentLine, line));
				lastCoords[0] = coords[0];
				lastCoords[1] = coords[1];
				break;
			}
			case PathIterator.SEG_CLOSE: {
				final Line2D.Double currentLine = new Line2D.Double(coords[0], coords[1], firstCoords[0],
						firstCoords[1]);
				if (currentLine.intersectsLine(line))
					intersections.add(getIntersection(currentLine, line));
				break;
			}
			default: {
				throw new Exception("Unsupported PathIterator segment type.");
			}
			}
			polyIt.next();
		}
		return intersections;

	}

	public static Point2D getIntersection(final Line2D.Double line1, final Line2D.Double line2) {

		final double x1, y1, x2, y2, x3, y3, x4, y4;
		x1 = line1.x1;
		y1 = line1.y1;
		x2 = line1.x2;
		y2 = line1.y2;
		x3 = line2.x1;
		y3 = line2.y1;
		x4 = line2.x2;
		y4 = line2.y2;
		final double x = ((x2 - x1) * (x3 * y4 - x4 * y3) - (x4 - x3) * (x1 * y2 - x2 * y1))
				/ ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));
		final double y = ((y3 - y4) * (x1 * y2 - x2 * y1) - (y1 - y2) * (x3 * y4 - x4 * y3))
				/ ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));

		return new Point2D.Double(x, y);

	}

	@Override
	public void run(String arg0) {
		// TODO Auto-generated method stub

	}

}
