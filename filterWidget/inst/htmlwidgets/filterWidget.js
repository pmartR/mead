HTMLWidgets.widget({

  name: 'filterWidget',

  type: 'output',

  factory: function(el, width, height) {

    // create an empty chart
    var chart = null;

    return {

      renderValue: function(x) {
          var toType = function(obj) {
            return ({}).toString.call(obj).match(/\s([a-zA-Z]+)/)[1].toLowerCase()
          }

            var arr = x.dataset[x.colName];
            var o_converter = d3.scale.ordinal();
            function binner(ars) {
                var a = [], b = [], prev;

                ars.sort();
                for ( var i = 0; i < ars.length; i++ ) {
                    if ( ars[i] !== prev ) {
                        a.push(ars[i]);
                        b.push(1);
                    } else {
                        b[b.length-1]++;
                    }
                    prev = ars[i];
                }

                return [a, b];
            }
            console.log(toType(arr[1]));
            console.log(binner(arr));
            console.log(x.dataset);
            console.log(x.colName);
            console.log(x.dataset[x.colName]);
            var bins = d3.layout.histogram()  // create layout object
                .bins(10)
                (arr);          // group the data into the bins
            var counts = [];

            for (i=0; i < bins.length; i++){
              counts[i] = bins[i].y;
            }
            if (toType(arr[1]) == "string"){
              counts = binner(arr)[1];
            }
            var bin_labs = [];
            for (i=0; i < bins.length; i++){
              bin_labs[i] = bins[i].x;
            }
            console.log(bin_labs);
            console.log(counts);
            var binDat = [];
            for (i=0; i < bins.length; i++){
              binDat.push({
                key: bin_labs[i],
                value: counts[i]
              });
            }
            console.log(binDat);

        // if the chart does not exist, create it via c3.generate
        if(chart===null){
            chart = c3.generate({
              // specify the container element we want the chart to render in
                bindto: el,
                data: {
                      // intialize with an empty array
                   columns: [
                     counts,
                   ],
                    // set chart types
                    type: 'bar',
                },
                axes: {
                  counts: bin_labs,
                },
                axis: {
                  y: {
                    label: x.colName,
                  }
                },
                // display a subchart - this will be used for brushing in a later stage
                subchart: {
                    show: true
                }
            });
            el.chart = chart;
        }
        el.chart.load({columns: ['y', counts]});



        // at this stage the chart always exists
        // get difference in keys
        //var old_keys = _.keys(chart.x());
        //var new_keys = _.keys(x.dataset);
        //var diff     = _.difference(old_keys,new_keys);

        // update the data and colors

      },
      resize: function(width, height) {

        // TODO: code to re-render the widget with a new size

      }

    };

  }
});
