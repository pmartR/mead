HTMLWidgets.widget({

  name: 'filterWidget',

  type: 'output',

  factory: function(el, width, height) {

    return {

      renderValue: function(x) {
        var arr = x.dataset[x.colName];

        var formatCount = d3.format(",.0f");

        var margin = {top: 10, right: 30, bottom: 30, left: 30},
        w = width - margin.left - margin.right,
        h = height - margin.top - margin.bottom;

        var svg = d3.select(el).append("svg")
          .style("width", width)
          .style("height", height);

        var g = svg.append("g")
          .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

        var x = d3.scaleLinear()
          .rangeRound([0, w]);

        var bins = d3.histogram()
            .domain(x.domain())
            .thresholds(x.ticks(10))
            (arr);

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
            .attr("height", function(d) { return h - y(d.length); });

        bar.append("text")
            .attr("dy", ".75em")
            .attr("y", 6)
            .attr("x", (x(bins[0].x1) - x(bins[0].x0)) / 2)
            .attr("text-anchor", "middle")
            .text(function(d) { return formatCount(d.length); });

        g.append("g")
            .attr("class", "axis axis--x")
            .attr("transform", "translate(0," + h + ")")
            .call(d3.axisBottom(x));


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
