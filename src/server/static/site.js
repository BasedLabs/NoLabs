var downloadPdbFile = function (proteinFileContent) {
    const obj = {
        bind: () => {
            $('#downloadPdbFile').on('click', () => {
                const filename = 'protein.pdb';
                const blob = new Blob([proteinFileContent], {type: 'text/plain'});
                if (window.navigator.msSaveOrOpenBlob) {
                    window.navigator.msSaveBlob(blob, filename);
                } else {
                    const elem = window.document.createElement('a');
                    elem.href = window.URL.createObjectURL(blob);
                    elem.download = filename;
                    document.body.appendChild(elem);
                    elem.click();
                    document.body.removeChild(elem);
                }
            });
        }
    }

    return obj;
}

var localisation = function (probabilities) {
    if (!document.getElementById('img'))
        throw Error('Declare img tag');
    const _createImage = (src) => {
        const image = new Image(500, 500);
        image.src = src;
        return image;
    };
    const img = document.getElementById('img');
    const obj = {
        hoverControlHighlight: {
            onmouseover: (el) => {
                const jqueryEl = $(el.target);
                jqueryEl.addClass('active');
                const localisationObjectKey = jqueryEl.data('localisationObjectKey');
                const highlightElement = obj.highlightData[localisationObjectKey];
                img.src = highlightElement.image.src;
            },
            onmouseleave: (el) => {
                $(el.target).removeClass('active');
                img.src = obj.highlightData.original.image.src;
            },
        },
        highlightData: {
            mithochondria: {
                image: _createImage('./static/images/mithochondria.png'),
                controlElementId: 'mithochondria-list-item',
                text: 'Mithochondria'
            },
            nucleus: {
                image: _createImage('./static/images/nucleus.png'),
                controlElementId: 'nucleus-list-item',
                text: 'Nucleus'
            },
            cytoplasm: {
                image: _createImage('./static/images/cytoplasm.png'),
                controlElementId: 'cytoplasm-list-item',
                text: 'Cytoplasm'
            },
            other: {
                image: _createImage('./static/images/original.png'),
                controlElementId: 'other-proteins-item',
                text: 'Other proteins'
            },
            extracellular: {
                image: _createImage('./static/images/original.png'),
                controlElementId: 'extracellular-proteins-item',
                text: 'Extracellular proteins'
            },
            original: {
                image: _createImage('./static/images/original.png')
            },
        },
        render: function () {
            const highlightElements = [];
            for (const key in obj.highlightData) {
                const data = obj.highlightData[key];
                const elId = data.controlElementId;
                const el = $('#' + elId);
                el.empty();
                if (probabilities[key] === undefined)
                    probabilities[key] = 0.000
                el.append(data.text + ' ' + (probabilities[key].toFixed(2) * 100) + '%');
                el.hover(obj.hoverControlHighlight.onmouseover,
                    obj.hoverControlHighlight.onmouseleave);
                el.data('localisationObjectKey', key);
                highlightElements.push(el);
            }

            img.src = obj.highlightData.original.image.src;
        }
    }
    return obj;
}

var solubility = function (solubility) {
    const obj = {
        render() {
            $('#solubilityText').text(`This protein is soluble with ${+(Math.round(solubility * 100.0 + "e+2") + "e-2")}% probability`);
        }
    }
    return obj;
}
var folding = function (proteinFileContent) {
    window.foldingInstance = undefined;
    $('#viewport').empty();
    $('#nav3dViewerTab').trigger('click');
    const stage = new NGL.Stage("viewport");
    stage.setParameters({backgroundColor: 'white'})
    const obj = {
        views: {
            default: 'default',
            cartoon: 'cartoon',
            backbone: 'backbone',
            ballsAndSticks: 'ball+stick',
            contact: 'contact',
            helixorient: 'helixorient',
            hyperball: 'hyperball',
            licorice: 'licorice',
            ribbon: 'ribbon',
            rope: 'rope',
            surface: 'surface',
            spacefill: 'spacefill',
            unitcell: 'unitcell'
        },
        reload: () => {
            stage.removeAllComponents();
            const proteinFileContentBlob = new Blob([proteinFileContent], {type: 'text/plain'});
            const proteinFile = new File([proteinFileContentBlob], 'protein.pdb', {type: 'text/plain'});
            stage.loadFile(proteinFile, {defaultRepresentation: true}).then((component) => {
                obj.component = component;
            });
        },
        setView: (viewName) => {
            if (viewName === obj.views.default) {
                obj.reload();
                return;
            }
            obj.component.removeAllRepresentations();
            obj.component.addRepresentation(viewName);
        },
        render: () => {
            obj.reload();
            window.foldingInstance = obj;
        }
    }
    return obj;
}

var drugTarget = function (drugTargetDiscoveryData) {
    window.drugTargetInstance = undefined;
    $('#viewport').empty();
    $('#nav3dViewerTab').trigger('click');
    const stage = new NGL.Stage("viewport");
    stage.setParameters({backgroundColor: 'white'})
    const obj = {
        currentPdb: undefined,
        currentSdf: undefined,
        views: {
            default: 'default',
            cartoon: 'cartoon',
            backbone: 'backbone',
            ballsAndSticks: 'ball+stick',
            contact: 'contact',
            helixorient: 'helixorient',
            hyperball: 'hyperball',
            licorice: 'licorice',
            ribbon: 'ribbon',
            rope: 'rope',
            surface: 'surface',
            spacefill: 'spacefill',
            unitcell: 'unitcell'
        },
        reload: () => {
            stage.removeAllComponents();
            const proteinFileContentBlob = new Blob([obj.currentPdb], {type: 'text/plain'});
            const ligandFileContentBlob = new Blob([obj.currentSdf], {type: 'text/plain'});
            const proteinFile = new File([proteinFileContentBlob], 'protein.pdb', {type: 'text/plain'});
            const sdfFile = new File([ligandFileContentBlob], 'ligand.sdf', {type: 'text/plain'});
            stage.loadFile(proteinFile, {defaultRepresentation: true}).then((component) => {
                obj.component = component;
            });

            stage.loadFile(sdfFile, {defaultRepresentation: true}).then((ligandComponent) => {
                obj.ligandComponent = ligandComponent;
            });
        },
        setView: (viewName) => {
            if (viewName === obj.views.default) {
                obj.reload();
                return;
            }
            obj.component.removeAllRepresentations();
            obj.component.addRepresentation(viewName);
            obj.ligandComponent.removeAllRepresentations();
            obj.ligandComponent.addRepresentation(viewName);
        },
        render: (pdb, sdf) => {
            if(!pdb || !sdf){
                pdb = drugTargetDiscoveryData[0].pdb;
                sdf = drugTargetDiscoveryData[0].sdf;
            }

            obj.currentPdb = pdb;
            obj.currentSdf = sdf;
            obj.reload();
            window.drugTargetInstance = obj;
        }
    }

    const tableData = [];
    for (let [index, val] of drugTargetDiscoveryData.entries()) {
        tableData.push({id: index, ligand: val.ligandName, protein: val.proteinName, affinity: val.affinity});
    }

    const rowSelect = (rowData) => {
        const drugTargetData = drugTargetDiscoveryData[rowData.id];
        obj.render(drugTargetData.pdb, drugTargetData.sdf);
    }

    $('#ligandProteinTable').bootstrapTable({
        data: tableData,
        onClickRow: (row, el, field) => {
            $('#ligandProteinTable tr').removeClass('active');
            $(el).addClass('active');
            rowSelect(row)
        }
    });

    $("#ligandProteinTableSearch").on("keyup", function () {
        const value = $(this).val().toLowerCase();
        $("#ligandProteinTable tr").filter(function () {
            $(this).toggle($(this).text().toLowerCase().indexOf(value) > -1)
        });
    });

    return obj;
}

var geneOntology = function (oboGraph) {
    const obj = {
        clear: () => {
            $('#geneOntologyContainer').empty();
        },
        render: () => {
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

                const graph = {name: "", children: Object.keys(oboGraph).map(key => convertToHierarchy(key, oboGraph))}
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

                node.append("text")
                    .on('mouseover', function (d, i) {
                        $(d.target).tooltip({
                            title: i.data.description
                        });
                        $(d.target).tooltip('show');
                    })
                    .on('mouseout', function (d, i) {
                        $(d.target).tooltip('hide');
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

    return obj;
}

var spinnerEnable = function () {
    $('#submitInference').attr('disabled', true);
    $('#spinner').removeClass('invisible');
}

var hideResultContainer = function () {
    $('#resultContainer').hide();
}

var showResultContainer = function () {
    $('#resultContainer').show();
}

