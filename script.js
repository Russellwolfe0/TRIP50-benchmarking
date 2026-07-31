function toggleCard(contentId, arrowId) {

    const content = document.getElementById(contentId);
    const arrow = document.getElementById(arrowId);

    content.classList.toggle("open");
    arrow.classList.toggle("rotated");

}
   
function maeColor(value, min, max) {

    let ratio = (value - min) / (max - min);

    if (max === min) {
        ratio = 0.5;
    }

    ratio = Math.max(0, Math.min(1, ratio));

    let r, g, b;

    if (ratio < 0.33) {

        // white -> pale green
        let t = ratio / 0.33;

        r = Math.round(255 - (25 * t));
        g = Math.round(255 - (5 * t));
        b = Math.round(255 - (35 * t));

    } 
    
    else if (ratio < 0.66) {

        // pale green -> pale yellow
        let t = (ratio - 0.33) / 0.33;

        r = Math.round(230 + (25 * t));
        g = Math.round(250 - (10 * t));
        b = Math.round(220 - (80 * t));

    } 
    
    else {

        // pale yellow -> pale red
        let t = (ratio - 0.66) / 0.34;

        r = Math.round(255);
        g = Math.round(240 - (70 * t));
        b = Math.round(140 + (80 * t));

    }

    return `rgb(${r}, ${g}, ${b})`;
}

document.addEventListener("DOMContentLoaded", function () {

    function toggleCard(contentId, arrowId) {

        const content = document.getElementById(contentId);
        const arrow = document.getElementById(arrowId);

        content.classList.toggle("open");
        arrow.classList.toggle("rotated");

    }


    Papa.parse("data/results.csv", {

        download: true,

        header: true,

        complete: function(results) {

            let thermoValues = results.data.map(
                row => parseFloat(row["MAE Thermo"])
            );
            
            let kineticValues = results.data.map(
                row => parseFloat(row["MAE Kinetics"])
            );
            
            let sumValues = results.data.map(
                row => parseFloat(row["Sum MAE"])
            );
            
            
            let thermoMin = Math.min(...thermoValues);
            let thermoMax = Math.max(...thermoValues);
            
            let kineticMin = Math.min(...kineticValues);
            let kineticMax = Math.max(...kineticValues);
            
            let sumMin = Math.min(...sumValues);
            let sumMax = Math.max(...sumValues);

            console.log(results.data);  // debugging

            let table = new DataTable('#results-table', {

                data: results.data,

                columns: [
                    { data: "Method" },
                    { data: "Type" },
                    {
                        data: "MAE Thermo",
                    
                        createdCell: function(cell, cellData) {
                    
                            let value = parseFloat(cellData);
                    
                            cell.style.backgroundColor = maeColor(
                                value,
                                0,
                                10
                            );
                    
                        }
                    },
                    {
                        data: "MAE Kinetics",
                    
                        createdCell: function(cell, cellData) {
                    
                            let value = parseFloat(cellData);
                    
                            cell.style.backgroundColor = maeColor(
                                value,
                                0,
                                10
                            );
                    
                        }
                    },
                    {
                        data: "Sum MAE",
                    
                        createdCell: function(cell, cellData) {
                    
                            let value = parseFloat(cellData);
                    
                            cell.style.backgroundColor = maeColor(
                                value,
                                0,
                                10
                            );
                    
                        }
                    },
                    {
                        data: "Average Time (s)",
                    
                        createdCell: function(cell, cellData) {
                    
                            let value = parseFloat(cellData);
                    
                            cell.style.backgroundColor = maeColor(
                                value,
                                0,
                                10
                            );
                    
                        }
                    },
                    {
                        data: "Efficiency Score",
                    
                        createdCell: function(cell, cellData) {
                    
                            let value = parseFloat(cellData);
                    
                            cell.style.backgroundColor = maeColor(
                                value,
                                0,
                                10
                            );
                    
                        }
                    },
                ]

            });

        }

    });

});