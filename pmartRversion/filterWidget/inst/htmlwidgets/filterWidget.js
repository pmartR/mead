HTMLWidgets.widget({

  name: 'filterWidget',

  type: 'output',

  factory: function(el, width, height) {
    var initialized = false;
    return {

      renderValue: function(val) {
        var arr = val.dataset[val.colName];
        if (!initialized) {
          initialized = true;
          if (HTMLWidgets.shinyMode) {
            Shiny.onInputChange("selected_data", arr);
            Shiny.onInputChange("brushed_on", val.colName);
          }
        }
        var ext = d3.extent(arr);

        var formatCount = d3.format(",.0f");

        var margin = {top: 10, right: 30, bottom: 30, left: 30},
        w = width - margin.left - margin.right,
        h = height - margin.top - margin.bottom;

        d3.select(el).append("p")
          .text(val.colName);
        var svg = d3.select(el).append("svg")
          .style("width", width)
          .style("height", height);

        var g = svg.append("g")
          .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

        var x = d3.scaleLinear()
          .rangeRound([0, w])
          .domain([ext[0],ext[1]].map(Math.round));

        var bins = d3.histogram()
            .domain(x.domain())
            // .thresholds(x.ticks(10))
            (arr);

        var rotate = d3.max(bins.map(function(d) {
          return (d.x0.toString()).length;
        })) > 1;

        var xaxis = d3.axisBottom(x)
          .ticks(bins.length);

        var y = d3.scaleLinear()
            .domain([0, d3.max(bins, function(d) { return d.length; })])
            .range([h, 0]);

        var bar = g.selectAll(".bar")
          .data(bins)
          .enter().append("g")
            .attr("class", "bar")
            .attr("transform", function(d) { return "translate(" + x(d.x0) + "," + y(d.length) + ")"; });

        bar.append("rect")
            .attr("x", 1)
            .attr("width", x(bins[0].x1) - x(bins[0].x0) - 1)
            // .attr("width", x(bins[0].x1) - x(bins[0].x0))
            .attr("height", function(d) { return h - y(d.length); });

        bar.append("text")
            .attr("dy", ".75em")
            .attr("y", 6)
            .attr("x", (x(bins[0].x1) - x(bins[0].x0)) / 2)
            .attr("text-anchor", "middle")
            .text(function(d) { return formatCount(d.length); });

        var bandwidth = bins[0].x1 - bins[0].x0;

        g.append("g")
            .attr("class", "axis axis--x")
            .attr("transform", "translate(0," + h + ")")
            .call(xaxis)
            .selectAll("text")
              .attr("y", function(d) { return rotate ? 5 : 10; })
              .attr("x", function(d) { return rotate ? 9 : 0; })
              .attr("transform", function(d) { return rotate ? "rotate(45)" : ""; })
              .style("text-anchor", function(d) { return rotate ? "start" : ""; });

      var brush = d3.brushX()
          .extent([[0, 0], [w, h]])
          .on("end", brushended);
          
      svg.append("g")
          .attr("class", "brush")
          .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
          .call(brush);
        function brushended() {
          if (!d3.event.sourceEvent) return; // Only transition after input.
          if (!d3.event.selection) return; // Ignore empty selections.
          var d0 = d3.event.selection.map(x.invert),
              d1 = d0.map(function(d) {
                var dd = Math.round(d * 100);
                var bw =  Math.round(bandwidth * 100);
                return (Math.round(dd/bw) * bw) / 100 ;
              });

          // d1.map(x.invert);
          // If empty or too small, transition to no brush
          if (d1[0] >= d1[1]) {
            d1[0] = d1[0];
            d1[1] = d1[0];
          }

          d3.select(this).transition().call(d3.event.target.move, d1.map(x));
          //return brushed data
          function checkSelection(arr) {
            return arr >= d1[0] && arr <= d1[1];
          }
          function isInRange(arr) {
            var indices = [], i;
            for(i = 0; i < arr.length; i++)
              if (arr[i] >= d1[0] && arr[i] <= d1[1])
                indices.push(i);
            return indices;
          }
          if (HTMLWidgets.shinyMode) {
              Shiny.onInputChange("selected_data", arr);
              Shiny.onInputChange("selected_data", arr.filter(checkSelection));
              Shiny.onInputChange("selected_indices", isInRange(arr));
              Shiny.onInputChange("brushed_on", val.colName);
          }

        }

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
