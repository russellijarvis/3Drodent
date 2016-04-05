<!DOCTYPE html>
<meta charset="utf-8">
<style>

body {
  font: 10px sans-serif;

}

.chord path {
  fill-opacity: .67;
  stroke: #000;
  stroke-width: .5px;
}

#tooltip{
    visibility: hidden;
    position: absolute;
    background-color: dodgerblue;
     border-radius: 4px;
     padding: 5px;
    z-index: 10;
    color:white;
    font-size:14px;
}
    

</style>


<body>  


<!--Offline  
<script type="text/javascript" src="d3.v3.min.js"></script> 

Online 
-->
<script src="//d3js.org/d3.v3.min.js"></script>
<script>

var width = 960,
    height = 500;

var color = d3.scale.category20();

var force = d3.layout.force()
    .charge(-120)
    .linkDistance(30)
    .size([width, height]);

var svg = d3.select("body").append("svg")
    .attr("width", width)
    .attr("height", height);



var hostname = window.location.hostname;
var isDev = (hostname == 'localhost') || (hostname == '127.0.0.1');
var tooltip = d3.select("#tooltip");//The Tooltip plugin is small pop-up box that appears when the user moves the mouse pointer over an element.

d3.json('js/global_whole_network.json', function(error, data){
  if (error) throw error;
  var matrix=data[2];
	  console.log(data);

  graph=data[0];
  //.push.apply(graph.nodes, graphi.nodes)
  graphi=data[1];//+data[1];
  //.push.apply(graph.nodes, graphi.nodes)
  entire_nodes=graph.nodes;
  entire_links=graph.links;
  
  entire_nodes.push.apply(graph.nodes, graphi.nodes);
  entire_links.push.apply(graph.links, graphi.links);

  var NameProvider = graph.nodes;
  force
      .nodes(entire_nodes)
      .links(entire_links)
      .start();

      //.nodes(graph.nodes)
      //.links(graph.links)

  var link = svg.selectAll(".link")
      .data(entire_links)
      //.data(graph.links)
    .enter().append("line")
      .attr("class", "link")
      .style("stroke-width", function(d) { return Math.sqrt(d.value); });

  var node = svg.selectAll(".node")
      .data(entire_nodes)
      //.data(graph.nodes)
    .enter().append("circle")
      .attr("class", "node")
      .attr("r", 5)
      .style("fill", function(d) { return color(d.group); })
      .call(force.drag);

  node.append("title")
      .text(function(d) { return d.name; });

  force.on("tick", function() {
    link.attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

    node.attr("cx", function(d) { return d.x; })
        .attr("cy", function(d) { return d.y; });


  });
});
/*Sums up to exactly 100*/

</script>
</div>