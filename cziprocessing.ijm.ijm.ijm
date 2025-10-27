// ==== CONFIGURATION ====

// Directory setup
inputDir = getDirectory("Choose folder with CZI files");
outputDir = getDirectory("Choose output folder");
list = getFileList(inputDir);

// === PREVIEW MODE ===
// If true → processes only the FIRST .czi file (for testing)
previewMode = false; // <-- set to false for full batch run

// Channel-specific rolling background values
rolling1 = 20; // DAPI (Blue)
rolling2 = 40; // Amyloid (Green)
rolling3 = 60; // IgG (Red)
rolling4 = 40; // AF594 (Purple)
rolling5 = 60; // Cy5 (Far Red)
saturated = 0.05; // Gentle contrast enhancement

// Bleed-through correction: how much Cy5 to subtract from AF594
bleedFraction = 0.10; // 15% of Cy5 subtracted from AF594 (adjust 0.1–0.25 as needed)

// ==== MAIN PROCESSING LOOP ====
for (i=0; i<list.length; i++) {
    if (endsWith(list[i], ".czi")) {

        // === Preview mode: process only first image ===
        if (previewMode && i > 0) break;

        print("Processing: " + list[i]);

        run("Bio-Formats Importer", "open=[" + inputDir + list[i] + "]");
        run("Split Channels");

        // Track available channels
        available = newArray(5);
        for (c=1; c<=5; c++) {
            title = "C" + c + "-" + list[i];
            if (isOpen(title)) available[c-1] = 1; else available[c-1] = 0;
        }

        // Assign titles
        c1Title = "C1-" + list[i];
        c2Title = "C2-" + list[i];
        c3Title = "C3-" + list[i];
        c4Title = "C4-" + list[i];
        c5Title = "C5-" + list[i];

        // === Background subtraction and contrast enhancement ===
        if (available[0]) { selectWindow(c1Title); run("Subtract Background...", "rolling=" + rolling1); run("Enhance Contrast", "saturated=" + saturated); }
        if (available[1]) { selectWindow(c2Title); run("Subtract Background...", "rolling=" + rolling2); run("Enhance Contrast", "saturated=" + saturated); }
        if (available[2]) { selectWindow(c3Title); run("Subtract Background...", "rolling=" + rolling3); run("Enhance Contrast", "saturated=" + saturated); }
        if (available[3]) { selectWindow(c4Title); run("Subtract Background...", "rolling=" + rolling4); run("Enhance Contrast", "saturated=" + saturated); }
        if (available[4]) { selectWindow(c5Title); run("Subtract Background...", "rolling=" + rolling5); run("Enhance Contrast", "saturated=" + saturated); }

        // === Bleed-through correction (AF594 - Cy5) ===
        if (available[3] && available[4]) {
            selectWindow(c4Title); run("Duplicate...", "title=AF594_temp");
            selectWindow(c5Title); run("Duplicate...", "title=Cy5_temp");

            selectWindow("Cy5_temp");
            run("Multiply...", "value=" + bleedFraction);

            selectWindow("AF594_temp");
            imageCalculator("subtract create", "AF594_temp", "Cy5_temp");
            rename("AF594_corrected");

            close(c4Title);
            c4Title = "AF594_corrected";
            close("AF594_temp");
            close("Cy5_temp");
        }

        // === Build merge command dynamically ===
        mergeCommand = "";
        colorList = newArray("B&W", "Green", "Red", "Violet", "Cyan");
        mergedColors = newArray();

        chanIndex = 1;
        for (c=1; c<=5; c++) {
            if (available[c-1]) {
                mergeCommand += "c" + chanIndex + "=";
                if (c==1) mergeCommand += c1Title + " ";
                else if (c==2) mergeCommand += c2Title + " ";
                else if (c==3) mergeCommand += c3Title + " ";
                else if (c==4) mergeCommand += c4Title + " ";
                else if (c==5) mergeCommand += c5Title + " ";
                mergedColors[chanIndex-1] = colorList[c-1];
                chanIndex++;
            }
        }
        mergeCommand += "create";
        run("Merge Channels...", mergeCommand);

        // === Apply display colors ===
        colorCommand = "channels=" + (chanIndex-1) + " colors=[";
        for (j=0; j<(chanIndex-1); j++) {
            if (j>0) colorCommand += ", ";
            colorCommand += mergedColors[j];
        }
        colorCommand += "] display=Composite";
        run("Channels...", colorCommand);

        // === Save composite ===
        mergedTitle = getTitle();
        selectWindow(mergedTitle);
        saveAs("Tiff", outputDir + "processed_" + replace(list[i], ".czi", "") + "_composite.tif");

        if (previewMode) {
            print("Preview mode: showing composite for inspection.");
            waitForUser("Preview complete. Check the image and click OK to continue or stop.");
        }

        run("Close All");
    }
}

print("Processing complete.");
if (previewMode) print("⚠️ Preview mode was ON — only the first file was processed.");
