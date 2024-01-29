<script lang="ts">
import {defineComponent, PropType} from 'vue'
import {AminoAcid} from "src/features/aminoAcid/geneOntology/types";

type RowType = {
  id: string;
  name: string;
  namespace: string;
}

export default defineComponent({
  name: "GeneOntologyTree",
  props: {
    oboGraph: {
      type: Object as PropType<{
        [id: string]: { name: string, namespace: string, edges: { [name: string]: Array<string> } }
      }>,
      required: true
    }
  },
  computed: {
    rows(): Array<RowType> {
      const rows = [];
      for (const id of Object.keys(this.oboGraph)) {
        rows.push({
          id: id,
          name: this.oboGraph[id].name,
          namespace: this.oboGraph[id].namespace,
        })
      }
      return rows;
    }
  },
  data() {
    return {
      filter: '',
      selectedRow: [] as Array<RowType>,
      columns: [
        {
          name: 'id',
          label: 'Id',
          align: 'left',
          sortable: 'false',
          field: (row: RowType) => row.id,
        },
        {
          name: 'name',
          label: 'Name',
          align: 'left',
          sortable: true,
          field: (row: RowType) => row.name,
        },
        {
          name: 'namespace',
          label: 'Namespace',
          align: 'left',
          sortable: true,
          field: (row: RowType) => row.namespace,
        }
      ]
    }
  },
  mounted() {
    this.render(this.oboGraph);
  },
  methods: {
    filterMethod(rows: Array<RowType>, terms: string): Array<RowType> {
      return rows.filter(x => x.name.toLowerCase().indexOf(terms) >= 0 ||
          x.namespace.toLowerCase().indexOf(terms) >= 0 ||
          x.id.toLowerCase().indexOf(terms) >= 0);
    },
    resetSelection() {
      d3.selectAll('circle').attr('r', 6);
    },
    geneOntologyTableRowClick(row: RowType) {
      this.resetSelection();
      this.selectedRow = [row];
      d3.select(`#circle-${row.id.replace('GO:', '')}`)
          .attr("r", 30);
    },
    render(oboGraph: {
      [name: string]: { name: string, namespace: string, edges: { [name: string]: Array<string> } }
    }) {
      const nodesAfterSecondLevel = [];

      function convertToHierarchy(id: string, graph: {
        [name: string]: { name: string, namespace: string, edges: { [name: string]: Array<string> } }
      }, relation: string) {
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

      window.onmousemove = function (e) {
        var x = e.clientX,
            y = e.clientY;
        let tooltip = document.getElementById("tooltip");
        tooltip.style.top = (y - 40) + 'px';
        tooltip.style.left = (x - 100) + 'px';
      };

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
              if(d.target.data.name === 'GO:0007275' || d.target.data.name === 'GO:0048731'){
                if(d.target.data.relation === '' || d.target.data.relation === null || d.target.data.relation === undefined){
                  debugger;
                }
              }
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
        }

        function hideTooltip() {
          var tooltip = document.getElementById("tooltip");
          setTimeout(() => {
            tooltip.style.display = "none";
          }, 100);
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

      chart();
    }
  }
})
</script>

<template>
  <div class="row q-ma-sm">
    <div class="col-md-4">
      <q-table
          title="Gene ontology"
          :rows="rows"
          :columns="columns"
          row-key="id"
          v-model:selected="selectedRow"
          :rows-per-page-options="[3, 5, 10]"
          :filter="filter"
          :filter-method="filterMethod"
      >
        <template v-slot:top-right>
          <q-input borderless dense debounce="300" v-model="filter" placeholder="Search">
            <template v-slot:append>
              <q-icon name="search"/>
            </template>
          </q-input>
        </template>
        <template v-slot:body="props">
          <q-tr :props="props" @click="geneOntologyTableRowClick(props.row)">
            <q-td
                auto-width
                v-for="col in props.cols"
                :key="col.name"
                :props="props"
                class="hover-finger"
            >
              {{ col.value }}
            </q-td>
          </q-tr>
        </template>
      </q-table>
      <svg id="geneOntologyLegend" height="300px">
        <defs>
          <!-- A marker to be used as an arrowhead -->
          <marker id="arrow" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="6" markerHeight="6"
                  orient="auto-start-reverse">
            <path d="M 0 0 L 10 5 L 0 10 z"/>
          </marker>
        </defs>
      </svg>
    </div>
    <div class="col-md-8">
      <div id="tooltip" style="position: fixed; display: none; color: black;"></div>
      <div id="geneOntologyContainer" style="background-color: white"></div>
    </div>
  </div>
</template>
