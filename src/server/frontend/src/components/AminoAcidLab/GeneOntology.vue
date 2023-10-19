<script>
export default {
    props: ['experiment'],
    mounted() {
        const oboGraph = this.experiment.data.oboGraph;
        const nodesAfterSecondLevel = [];

        function convertToHierarchy(id, graph, relation) {
            const links = graph[id].edges;
            const name = graph[id].name;
            const namespace = graph[id].namespace;

            const children = [];

            if (links) {
                for (let link in links) {
                    children.push(convertToHierarchy(link, graph, links[link][0]));
                    const child = children[children.length - 1];
                    nodesAfterSecondLevel.push(child.name);
                }
            }

            return {
                name: id,
                description: name,
                namespace: namespace,
                relation: relation,
                children: children.length ? children : []
            };
        }

        const table = () => {
            const tableData = Object.keys(oboGraph).map(key => {
                return {
                    id: key,
                    name: oboGraph[key].name,
                    namespace: oboGraph[key].namespace
                }
            });

            const resetSelection = () => {
                d3.selectAll('circle').attr('r', 6);
            }

            const rowSelect = (rowData) => {
                resetSelection();
                d3.select(`#circle-${rowData.id.replace('GO:', '')}`)
                    .attr("r", 20);
            }

            $('#geneOntologyTable').bootstrapTable({
                data: tableData,
                onClickRow: (row, el, field) => {
                    $('#geneOntologyTable tr').removeClass('active');
                    $(el).addClass('active');
                    rowSelect(row)
                }
            });

            $("#geneOntologyTableSearch").on("keyup", function () {
                resetSelection();
                const value = $(this).val().toLowerCase();
                $("#geneOntologyTable tr").filter(function () {
                    $(this).toggle($(this).text().toLowerCase().indexOf(value) > -1)
                });
            });
        }

        const chart = () => {
            const width = 1024;

            const graph = { name: "", children: Object.keys(oboGraph).map(key => convertToHierarchy(key, oboGraph)) }
            graph.children = graph.children.filter(x => !nodesAfterSecondLevel.includes(x.name));

            // Compute the tree height; this approach will allow the height of the
            // SVG to scale according to the breadth (width) of the tree layout.
            const root = d3.hierarchy(graph);
            const dx = 40;
            const dy = width / (root.height + 1);

            // Create a tree layout.
            const tree = d3.tree().nodeSize([dx, dy]);

            // Sort the tree and apply the layout.
            root.sort((a, b) => d3.ascending(a.data.name, b.data.name));
            tree(root);

            // Compute the extent of the tree. Note that x and y are swapped here
            // because in the tree layout, x is the breadth, but when displayed, the
            // tree extends right rather than down.
            let x0 = Infinity;
            let x1 = -x0;
            root.each(d => {
                if (d.x > x1) x1 = d.x;
                if (d.x < x0) x0 = d.x;
            });

            // Compute the adjusted height of the tree.
            const height = x1 - x0 + dx * 2;

            const svg = d3.create("svg")
                .attr("width", width)
                .attr("height", height)
                .attr("viewBox", [-dy / 3, x0 - dx, width, height])
                .attr("style", "max-width: 100%; height: auto; font: 15px sans-serif;");

            const pathIdFactory = (d) => {
                return d.source.data.name.replace('GO:', '') + '-' + d.target.data.name.replace('GO:', '');
            }

            const circleIdFactory = (d) => {
                return 'circle-' + d.data.name.replace('GO:', '');
            }

            const textIdFactory = (d) => {
                return 'text-' + d.data.name.replace('GO:', '');
            }

            const namespaceColorFactory = (d) => {
                if (typeof d === 'string') {
                    if (d === 'Biological process')
                        return '#1ddb05';

                    if (d === 'Cellular component')
                        return '#244fff';

                    if (d === 'Molecular function')
                        return '#ff0f0f';
                }

                if (d.data.namespace === 'biological_process')
                    return '#1ddb05';

                if (d.data.namespace === 'cellular_component')
                    return '#244fff';

                if (d.data.namespace === 'molecular_function')
                    return '#ff0f0f';
            }

            const link = svg.append('g')
                .selectAll()
                .data(root.links())
                .join("path")
                .attr("d", d3.linkHorizontal()
                    .x(d => d.y)
                    .y(d => d.x))
                .attr('id', pathIdFactory)
                .each(function (d, i) {
                    d3.select(this)
                        .attr("fill", "none")
                        .attr("stroke", "#555")
                        .attr("stroke-opacity", 0.4)
                        .attr("stroke-width", 3);
                });

            svg.append('g').selectAll()
                .data(root.links())
                .join('text')
                .attr('font-family', 'Verdana')
                .attr('font-size', 14)
                .attr('fill', 'black')
                .append('textPath')
                .each(function (d, i) {
                    const href = pathIdFactory(d);
                    d3.select(this)
                        .attr('xlink:xlink:href', '#' + href)
                        .attr('text-anchor', 'middle')
                        .attr('startOffset', '50%')
                        .text(d.target.data.relation)
                });

            const node = svg.append("g")
                .attr("stroke-linejoin", "round")
                .attr("stroke-width", 6)
                .selectAll()
                .data(root.descendants())
                .join("g")
                .attr("transform", d => `translate(${d.y},${d.x})`);

            node.append("circle")
                .attr("fill", d => {
                    return namespaceColorFactory(d);
                })
                .attr("r", 6)
                .attr('id', circleIdFactory);

            function showTooltip(evt, text) {
                let tooltip = document.getElementById("tooltip");
                tooltip.innerHTML = text;
                tooltip.style.display = "block";
                tooltip.style.left = evt.pageX + 10 + 'px';
                tooltip.style.top = evt.pageY + 10 + 'px';
            }

            function hideTooltip() {
                var tooltip = document.getElementById("tooltip");
                tooltip.style.display = "none";
            }

            node.append("text")
                .on('mouseover', function (d, i) {
                    showTooltip(d, i.data.description);
                })
                .on('mouseout', function (d, i) {
                    hideTooltip();
                }).attr("dy", "0.40em")
                .attr("x", d => d.children ? -6 : 6)
                .attr("text-anchor", d => d.children ? "end" : "start")
                .text(d => d.data.name)
                .clone(true).lower()
                .attr("stroke", "white")
                .attr('id', textIdFactory);

            var legend = d3.select("#geneOntologyLegend")
            var keys = ["Biological process", "Cellular component", "Molecular function"]
            legend.selectAll("mydots")
                .data(keys)
                .enter()
                .append("circle")
                .attr("cx", 100)
                .attr("cy", function (d, i) {
                    return 100 + i * 25
                }) // 100 is where the first dot appears. 25 is the distance between dots
                .attr("r", 7)
                .style("fill", function (d) {
                    return namespaceColorFactory(d)
                })

            legend.selectAll("mylabels")
                .data(keys)
                .enter()
                .append("text")
                .attr("x", 120)
                .attr("y", function (d, i) {
                    return 100 + i * 25
                })
                .style("fill", function (d) {
                    return namespaceColorFactory(d)
                })
                .text(function (d) {
                    return d
                })
                .attr("text-anchor", "left")
                .style("alignment-baseline", "middle")

            $('#geneOntologyContainer').append(svg.node());
        };

        table();
        chart();
    }
}
</script>

<template>
    <div class="col-md-12 row" id="geneOntology">
        <div class="row">
            <div class="col-md-2">
                <svg id="geneOntologyLegend" height="300px">
                    <defs>
                        <!-- A marker to be used as an arrowhead -->
                        <marker id="arrow" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="6" markerHeight="6"
                            orient="auto-start-reverse">
                            <path d="M 0 0 L 10 5 L 0 10 z" />
                        </marker>
                    </defs>
                </svg>
            </div>
            <div id="geneOntologyContainer" class="col-md-10"></div>
        </div>
        <div class="row">
            <div class="col-md-2"></div>
            <div class="col-md-8 container">
                <div id="tooltip" display="none" style="position: absolute; display: none;"></div>
                <label for="geneOntologyTableSearch">Search for Gene ID or description, click to
                    highlight:</label><input class="form-control" id="geneOntologyTableSearch" type="text"
                    placeholder="Search..">
                <br>
                <table style="max-height: 200px; overflow-y:scroll;" class="table table-bordered table-striped"
                    id="geneOntologyTable">
                    <thead>
                        <tr>
                            <th data-field="id">Gene ID</th>
                            <th data-field="name">Name</th>
                            <th data-field="namespace">Namespace</th>
                        </tr>
                    </thead>
                </table>
            </div>
            <div class="col-md-2"></div>
        </div>
    </div>
</template>
