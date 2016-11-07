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

            var bins = d3.layout.histogram()  // create layout object
                .bins(10)
                (arr);          // group the data into the bins
            var counts = [];

            for (i=0; i < bins.length; i++){
              counts[i] = bins[i].y;
            }

            var bin_labs = [];
            for (i=0; i < bins.length; i++){
              bin_labs[i] = Math.round(bins[i].x*1000)/1000;
            }

            if (toType(arr[1]) == "string"){
              tallies = binner(arr);
              counts = tallies[1];
              bin_labs = tallies[0];
            }

            var binDat = [];
            for (i=0; i < bin_labs.length; i++){
              binDat.push({
                bin_labels: bin_labs[i],
                bin_counts: counts[i]
              });
            }


        // if the chart does not exist, create it via c3.generate
        if(chart===null){
            chart = c3.generate({
              // specify the container element we want the chart to render in
                bindto: el,

                data: {
                      // intialize with an empty array
                  json: binDat,
                  keys: {
                    x: 'bin_labels',
                    value: ['bin_counts'],
                  },
                    // set chart types
                  type: 'bar',

                },
                 legend: {
                    show: false
                  },
                axis: {
                  x: {
                    type: 'category',
                  },
                  y: {
                    label: {
                      text: x.colName,
                      position: 'outer-middle'
                    }
                  }
                },
                // display a subchart - this will be used for brushing in a later stage
                subchart: {
                    show: true,
                    //onbrush: function (domain) {
                      //console.log(this, domain);
                    //}
                }
            });
            el.chart = chart;
        }
        el.chart.load(binDat);



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
