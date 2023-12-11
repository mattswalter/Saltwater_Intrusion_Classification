///// Classification of land cover in Delmarva based on NAIP and Landsat data
///// Author: Matthew Walter, mswalter@udel.edu

///// Inputs needed: Date range for NAIP and Landsat images, 
///// shapefile of study area (studyarea), 
///// point file with reference data for classifier (ref_points), 
///// geometry covering entire study area for export (geometry)

// Import NAIP data and filter by date and study area
var dataset = ee.ImageCollection('USDA/NAIP/DOQQ') 
.filterBounds(studyarea) 
.filterDate('2017-01-01','2017-09-30')
.max()
.clip(studyarea);
var Newclip = dataset.clip(studyarea);
Map.addLayer(dataset,{},'NAIP data');


// Define a boxcar kernel for smoothing
var boxcar = ee.Kernel.square({
  radius: 3, units: 'pixels', magnitude: 1
});

// Smooth the image by convolving with the boxcar kernel.
var smooth = dataset.convolve(boxcar);
Map.addLayer(smooth, {min: 100, max: 200}, 'NAIP smoothed');
var smooth_clip = smooth.clip(studyarea)

// Calculate NDVI from NAIP
var nir = Newclip.select('N');
var red = Newclip.select('R');
var green = Newclip.select('G');
var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');
var dataset = dataset.addBands(ndvi)

// Calculate additional PC bands from NAIP and NDVI
var image = dataset
// Display the input imagery and the region in which to do the PCA.
var region = image.geometry();

// Set some information about the input to be used later.
var scale = 30;
var bandNames = image.bandNames();

// Mean center the data to enable a faster covariance reducer
// and an SD stretch of the principal components.
var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: studyarea,
    scale: scale,
    maxPixels: 1e9
});
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = image.subtract(means);

// This helper function returns a list of new band names.
var getNewBandNames = function(prefix) {
  var seq = ee.List.sequence(1, bandNames.length());
  return seq.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).int());
  });
};

// This function accepts mean centered imagery, a scale and
// a region in which to perform the analysis.  It returns the
// Principal Components (PC) in the region as a new image.
var getPrincipalComponents = function(centered, scale, region) {
  // Collapse the bands of the image into a 1D array per pixel.
  var arrays = centered.toArray();

  // Compute the covariance of the bands within the region.
  var covar = arrays.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: studyarea,
    scale: scale,
    maxPixels: 1e9
  });

  // Get the 'array' covariance result and cast to an array.
  // This represents the band-to-band covariance within the region.
  var covarArray = ee.Array(covar.get('array'));

  // Perform an eigen analysis and slice apart the values and vectors.
  var eigens = covarArray.eigen();

  // This is a P-length vector of Eigenvalues.
  var eigenValues = eigens.slice(1, 0, 1);
  // This is a PxP matrix with eigenvectors in rows.
  var eigenVectors = eigens.slice(1, 1);

  // Convert the array image to 2D arrays for matrix computations.
  var arrayImage = arrays.toArray(1);

  // Left multiply the image array by the matrix of eigenvectors.
  var principalComponents = ee.Image(eigenVectors).matrixMultiply(arrayImage);

  // Turn the square roots of the Eigenvalues into a P-band image.
  var sdImage = ee.Image(eigenValues.sqrt())
    .arrayProject([0]).arrayFlatten([getNewBandNames('sd')]);

  // Turn the PCs into a P-band image, normalized by SD.
  return principalComponents
    // Throw out an an unneeded dimension, [[]] -> [].
    .arrayProject([0])
    // Make the one band array image a multi-band image, [] -> image.
    .arrayFlatten([getNewBandNames('pc')])
    // Normalize the PCs by their SDs.
    .divide(sdImage);
};

// Get the PCs at the specified scale and in the specified region
var pcImage = getPrincipalComponents(centered, scale, region);

for (var i = 0; i < bandNames.length().getInfo(); i++) {
  var band = pcImage.bandNames().get(i).getInfo();}
print(pcImage)

var pc3 = pcImage.select(['pc3']).rename('PC3')

var pc1 = pcImage.select(['pc1']).rename('PC1')

var pc2 = pcImage.select(['pc2']).rename('PC2')

var pc4 = pcImage.select(['pc4']).rename('PC4')

var pc5 = pcImage.select(['pc5']).rename('PC5')

var pca_1 = pc1.addBands(pc2).addBands(pc3).addBands(pc4).addBands(pc5)

// Calculate additional indices from NAIP: ndvi, dvi, ndwi, si, rvi
var nir = Newclip.select('N');
var red = Newclip.select('R');
var green = Newclip.select('G');
var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');
var dvi = nir.subtract(red).rename('DVI');
var rvi = red.divide(nir).rename('RVI');
var ndwi = green.subtract(nir).divide(green.add(nir)).rename('NDWI');

var ndwi = Newclip.expression(
  '((Green - NIR)/(Green + NIR))',{
    'Green': Newclip.select('G'),
    'NIR': Newclip.select('N')
  }).rename('NDWI');
var si = Newclip.expression(
    '(256-BLUE)*(256-GREEN)', {
      'NIR': Newclip.select('N'),
      'RED': Newclip.select('R'),
      'BLUE': Newclip.select('B'),
      'GREEN': Newclip.select('G')
}).rename('SI');

// Combine all NAIP and NAIP-derived bands
var Newclip = Newclip.addBands(ndvi).addBands(si).addBands(ndwi).addBands(rvi).addBands(dvi).addBands(pc1).addBands(pc2).addBands(pc3).addBands(pc4).addBands(pc5).addBands(smooth_clip);


// Import Landsat 8 Imagery

// Cloud mask
var maskL8 = function(image) {
  var qa = image.select('BQA');
  /// Check that the cloud bit is off.
  // See https://landsat.usgs.gov/collectionqualityband
  var mask = qa.bitwiseAnd(1 << 4).eq(0);
  return image.updateMask(mask);
}

// Filter TOA image and calculate EVI for each season
var land = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA')
    .filterDate('2017-01-01', '2017-12-31')
    .filterBounds(studyarea)
    .map(maskL8)
    .select('B2', 'B4', 'B5')

// Calculate season 1 and smooth    
var s21 = ee.Image(land.filterDate('2017-04-01', '2017-06-30').median().clip(studyarea));
var evi1 = s21.expression(
  '2.5 * ((NIR-RED) / (NIR + 6 * RED - 7.5* BLUE +1))', {
    NIR:s21.select('B5'),
    RED:s21.select('B4'),
    BLUE:s21.select('B2')
  }).rename('evi_spring');
var evi1 = evi1.convolve(boxcar);

// Calculate season 2 and smooth
var s22 = ee.Image(land.filterDate('2017-07-01', '2017-09-30').median().clip(studyarea));
var evi2 = s22.expression(
  '2.5 * ((NIR-RED) / (NIR + 6 * RED - 7.5* BLUE +1))', {
    NIR:s22.select('B5'),
    RED:s22.select('B4'),
    BLUE:s22.select('B2')
  }).rename('evi_summer');
var evi2 = evi2.convolve(boxcar)

// Calculate season 3 and smooth
var s23 = ee.Image(land.filterDate('2017-10-01', '2017-12-31').median().clip(studyarea));
var evi3 = s23.expression(
  '2.5 * ((NIR-RED) / (NIR + 6 * RED - 7.5* BLUE +1))', {
    NIR:s23.select('B5'),
    RED:s23.select('B4'),
    BLUE:s23.select('B2')
  }).rename('evi_fall');
var evi3 = evi3.convolve(boxcar)

// Calculate season 4 and smooth
var s24 = ee.Image(land.filterDate('2017-01-01', '2017-03-31').median().clip(studyarea));
var evi4 = s24.expression(
  '2.5 * ((NIR-RED) / (NIR + 6 * RED - 7.5* BLUE +1))', {
    NIR:s24.select('B5'),
    RED:s24.select('B4'),
    BLUE:s24.select('B2')
  }).rename('evi_winter');
var evi4 = evi4.convolve(boxcar);

// Import thermal bands from Landsat
var land = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA')
    .filterDate('2017-01-01', '2017-12-31')
    .filterBounds(studyarea)
    .map(maskL8)
    .select('B11', 'B10')

// Calculate the average thermal values for each season and smooth
var l1 = ee.Image(land.filterDate('2017-04-01', '2017-06-30').median().clip(studyarea));
var l1 = l1.convolve(boxcar);
var l2 = ee.Image(land.filterDate('2017-07-01', '2017-09-30').median().clip(studyarea));
var l2 = l2.convolve(boxcar);
var l3 = ee.Image(land.filterDate('2017-10-01', '2017-12-31').median().clip(studyarea))
var l3 = l3.convolve(boxcar);
var l4 = ee.Image(land.filterDate('2017-01-01', '2017-03-31').median().clip(studyarea));
var l4 = l4.convolve(boxcar);

// Merge all NAIP bands and Landsat bands
var Newclip = Newclip.addBands(evi1).addBands(evi2).addBands(evi3).addBands(evi4).addBands(l1).addBands(l2).addBands(l3).addBands(l4)

// Set field with land cover values from reference data
var label = 'Class'; 


// Randomly split all training points into testing (70%) and training (30%) for each class
var c1 = (ref_points.filter(ee.Filter.eq('Class', 1)));
var c1 = c1.randomColumn();
var train_c1 = c1.filter(ee.Filter.lt('random', 0.7));
var val_c1 = c1.filter(ee.Filter.gte('random', 0.7));

var c2 = (ref_points.filter(ee.Filter.eq('Class', 2)));
var c2 = c2.randomColumn();
var train_c2 = c2.filter(ee.Filter.lt('random', 0.7));
var val_c2 = c2.filter(ee.Filter.gte('random', 0.7));

var c3 = (ref_points.filter(ee.Filter.eq('Class', 3)));
var c3 = c3.randomColumn();
var train_c3 = c3.filter(ee.Filter.lt('random', 0.7));
var val_c3 = c3.filter(ee.Filter.gte('random', 0.7));

var c4 = (ref_points.filter(ee.Filter.eq('Class', 4)));
var c4 = c4.randomColumn();
var train_c4 = c4.filter(ee.Filter.lt('random', 0.7));
var val_c4 = c4.filter(ee.Filter.gte('random', 0.7));

var c5 = (ref_points.filter(ee.Filter.eq('Class', 5)));
var c5 = c5.randomColumn();
var train_c5 = c5.filter(ee.Filter.lt('random', 0.7));
var val_c5 = c5.filter(ee.Filter.gte('random', 0.7));

var c6 = (ref_points.filter(ee.Filter.eq('Class', 6)));
var c6 = c6.randomColumn();
var train_c6 = c6.filter(ee.Filter.lt('random', 0.7));
var val_c6 = c6.filter(ee.Filter.gte('random', 0.7));

var c7 = (ref_points.filter(ee.Filter.eq('Class', 7)));
var c7 = c7.randomColumn();
var train_c7 = c7.filter(ee.Filter.lt('random', 0.7));
var val_c7 = c7.filter(ee.Filter.gte('random', 0.7));

var c8 = (ref_points.filter(ee.Filter.eq('Class', 8)));
var c8 = c8.randomColumn();
var train_c8 = c8.filter(ee.Filter.lt('random', 0.7));
var val_c8 = c8.filter(ee.Filter.gte('random', 0.7));


// Merge training points
var all_merged_points_train = train_c1.merge(train_c2).merge(train_c3).merge(train_c4).merge(train_c5).merge(train_c6).merge(train_c7).merge(train_c8)
var training = Newclip.sampleRegions({
collection: all_merged_points_train,
properties: [label],
scale: 1 
});

// Train classifier and apply to dataset
var trained = ee.Classifier.smileRandomForest(100).train(training,label);
var classified = Newclip.classify(trained);

// Map classification
Map.addLayer(classified,
{min: 1, max: 8, palette: ['green','teal','purple','red','blue','yellow','tan','pink']},
'classification');

// Get a confusion matrix representing resubstitution accuracy.
var all_merged_points_val = val_c1.merge(val_c2).merge(val_c3).merge(val_c4).merge(val_c5).merge(val_c6).merge(val_c7).merge(val_c8)
var trainAccuracy = trained.confusionMatrix();
print('Resubstitution error matrix: ', trainAccuracy);
print('Training overall accuracy: ', trainAccuracy.accuracy());

// Sample the input with a different random seed to get validation data.
var validation = Newclip.sampleRegions({
    collection: all_merged_points_val,
    properties: [label],
    scale: 1  //scale is pixel
});

// Classify the validation data.
var validated = validation.classify(trained);

// Get a confusion matrix representing expected accuracy.
var testAccuracy = validated.errorMatrix('Class', 'classification');
print('Validation error matrix: ', testAccuracy);
print('Validation overall accuracy: ', testAccuracy.accuracy());

// Calculate and plot variable importance
var dict = trained.explain();
print('Explain:',dict);
 
var variable_importance = ee.Feature(null, ee.Dictionary(dict).get('importance'));
 
var chart =
ui.Chart.feature.byProperty(variable_importance)
.setChartType('ColumnChart')
.setOptions({
title: 'Random Forest Variable Importance',
legend: {position: 'none'},
hAxis: {title: 'Bands'},
vAxis: {title: 'Importance'}
});
 
print(chart);

// Export classified image and testing points
Export.image.toDrive({
  image: classified,
  description: 'swi_de_17_10_14',
  scale: 1,
  region: geometry,
  maxPixels: 10000000000000 
});

Export.table.toDrive({
  collection: all_merged_points_val,
  description: 'de_swi_points_11_9_20',
  fileFormat: 'SHP'
});
