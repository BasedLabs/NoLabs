<template>
  <q-separator></q-separator>
  <q-card>
    <div class="row q-gutter-md">
      <div class="col-12">
        <q-separator></q-separator>
        <q-card-section>
          <div class="text-h6">Job Details</div>
        </q-card-section>
        <q-card-section>
          <q-list>
            <q-item>
              <q-item-section>
                <q-item-label>Job ID</q-item-label>
              </q-item-section>
              <q-item-section>{{ jobId }}</q-item-section>
            </q-item>
            <q-item>
              <q-item-section>
                <q-item-label>Job Name</q-item-label>
              </q-item-section>
              <q-item-section>
                <q-input v-model="editableJobName" @keyup.enter="updateJobName" dense clearable />
                <q-btn v-if="editableJobName !== (job?.job_name || '')" icon="check" color="info" flat
                       @click="updateJobName" />
              </q-item-section>
            </q-item>
            <q-item>
              <q-item-section>
                <q-item-label>Job Status</q-item-label>
              </q-item-section>
              <q-item-section>
                <template v-if="jobStatus === null">
                  Loading...
                </template>
                <template v-else>
                  <q-spinner v-if="jobStatus.running" color="primary" size="20px" />
                  {{ jobStatusText }}
                </template>
              </q-item-section>
            </q-item>
          </q-list>
        </q-card-section>
      </div>
      <q-separator></q-separator>
      <div class="col-12">
        <q-separator></q-separator>
        <q-card-section>
          <div class="text-h6">Job Inputs</div>
        </q-card-section>
        <q-card-section>
          <q-list>
            <q-item>
              <q-item-section>
                <q-item-label>Protein FASTA Content</q-item-label>
              </q-item-section>
              <q-item-section class="fasta-content-container">
                <div class="fasta-content">{{ protein?.fasta_content }}</div>
              </q-item-section>
            </q-item>
            <q-item>
              <q-item-section>
                <q-item-label>Descriptions</q-item-label>
              </q-item-section>
              <q-item-section>
                <q-input v-model.number="jobInputs.descriptions" @change="updateJobInputs" dense type="number" />
              </q-item-section>
            </q-item>
            <q-item>
              <q-item-section>
                <q-item-label>Alignments</q-item-label>
              </q-item-section>
              <q-item-section>
                <q-input v-model.number="jobInputs.alignments" @change="updateJobInputs" dense type="number" />
              </q-item-section>
            </q-item>
            <q-item>
              <q-item-section>
                <q-item-label>Hitlist Size</q-item-label>
              </q-item-section>
              <q-item-section>
                <q-input v-model.number="jobInputs.hitlist_size" @change="updateJobInputs" dense type="number" />
              </q-item-section>
            </q-item>
            <q-item>
              <q-item-section>
                <q-item-label>Expect</q-item-label>
              </q-item-section>
              <q-item-section>
                <q-input v-model.number="jobInputs.expect" @change="updateJobInputs" dense type="number" />
              </q-item-section>
            </q-item>
          </q-list>
        </q-card-section>
      </div>
      <q-separator></q-separator>
      <div class="col-12">
        <q-separator></q-separator>
        <q-card-section>
          <div class="text-h6">Result</div>
        </q-card-section>
        <q-card-section>
          <!-- Query Info -->
          <div v-if="jobHasGeneratedData">
            <div v-for="result in job?.result" :key="result.query_id">
              <p>Program: {{ result.program }}</p>
              <p>Database: {{ result.database }}</p>
              <p>Query ID: {{ result.query_id }}</p>
              <p>Query Definition: {{ result.query_def }}</p>
              <p>Query Length: {{ result.query_len }}</p>
            </div>
            <!-- Hit Visualization -->
            <div id="hit-visualization" class="hit-visualization full-width"></div>
            <q-btn @click="downloadCSV" color="primary" label="Download Results as CSV" class="q-my-md" />
            <!-- Table with HSPs including Hit Info -->
            <q-table
              :rows="visibleHits"
              :columns="columns"
              row-key="id"
              class="q-mt-md"
            >
              <template v-slot:body-cell="props">
                <q-td :props="props" :style="{ maxWidth: '300px', overflowX: 'auto', whiteSpace: 'nowrap' }">
                  <div v-if="props.col.field === 'definition' || props.col.field === 'qseq' || props.col.field === 'hseq' || props.col.field === 'midline'" class="cell-content">
                    <div class="scrollable-cell-content">{{ props.row[props.col.field] }}</div>
                    <q-btn
                      flat
                      dense
                      round
                      icon="content_copy"
                      @click="copyToClipboard(props.row[props.col.field])"
                      v-ripple
                      class="copy-btn"
                      size="sm"
                    />
                  </div>
                  <div v-else>
                    {{ props.row[props.col.field] }}
                  </div>
                </q-td>
              </template>
            </q-table>
          </div>
          <div v-else>
            No results available.
          </div>
        </q-card-section>
      </div>
    </div>
  </q-card>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import { QSpinner, QInput, QBtn, QTable, QTooltip, QIcon, QTd } from 'quasar';
import * as d3 from 'd3';
import {
  nolabs__application__use_cases__blast__api_models__JobResponse,
  ProteinContentResponse,
  nolabs__application__use_cases__blast__api_models__GetJobStatusResponse,
  nolabs__application__use_cases__blast__api_models__SetupJobRequest
} from "../../../../refinedApi/client";
import {
  getBlastJobApi,
  getProteinContent,
  getJobStatus,
  changeJobName,
  setupBlastJob
} from "../../refinedApi";

export default defineComponent({
  name: 'BlastJob',
  props: {
    jobId: String,
  },
  data() {
    return {
      experimentId: null as string | null,
      job: null as nolabs__application__use_cases__blast__api_models__JobResponse | null,
      protein: null as ProteinContentResponse | null,
      jobStatus: null as nolabs__application__use_cases__blast__api_models__GetJobStatusResponse | null,
      editableJobName: '' as string,
      visibleHitId: null as string | null, // Store the currently visible hit ID
      jobInputs: {
        descriptions: null,
        alignments: null,
        hitlist_size: null,
        expect: null
      }
    };
  },
  computed: {
    jobHasGeneratedData(): boolean | null {
      return this.job && this.job.result.length > 0;
    },
    jobStatusText(): string {
      if (this.jobStatus === null) {
        return '';
      }
      return this.jobStatus.running ? 'Running...' : 'Not running';
    },
    columns() {
      return [
        { name: 'hit_id', label: 'Hit ID', field: 'hit_id', align: 'left' },
        { name: 'definition', label: 'Definition', field: 'definition', align: 'left', sortable: true },
        { name: 'accession', label: 'Accession', field: 'accession', align: 'left' },
        { name: 'length', label: 'Length', field: 'length', align: 'left' },
        { name: 'num', label: 'HSP Number', field: 'num', align: 'left' },
        { name: 'bit_score', label: 'Bit Score', field: 'bit_score', align: 'left' },
        { name: 'score', label: 'Score', field: 'score', align: 'left' },
        { name: 'evalue', label: 'E-value', field: 'evalue', align: 'left' },
        { name: 'query_from', label: 'Query From', field: 'query_from', align: 'left' },
        { name: 'query_to', label: 'Query To', field: 'query_to', align: 'left' },
        { name: 'hit_from', label: 'Hit From', field: 'hit_from', align: 'left' },
        { name: 'hit_to', label: 'Hit To', field: 'hit_to', align: 'left' },
        { name: 'query_frame', label: 'Query Frame', field: 'query_frame', align: 'left' },
        { name: 'hit_frame', label: 'Hit Frame', field: 'hit_frame', align: 'left' },
        { name: 'identity', label: 'Identity', field: 'identity', align: 'left' },
        { name: 'positive', label: 'Positive', field: 'positive', align: 'left' },
        { name: 'gaps', label: 'Gaps', field: 'gaps', align: 'left' },
        { name: 'align_len', label: 'Alignment Length', field: 'align_len', align: 'left' },
        { name: 'qseq', label: 'Query Sequence', field: 'qseq', align: 'left', sortable: true },
        { name: 'hseq', label: 'Hit Sequence', field: 'hseq', align: 'left', sortable: true },
        { name: 'midline', label: 'Midline', field: 'midline', align: 'left', sortable: true }
      ];
    },
    visibleHits() {
      if (this.visibleHitId) {
        return this.job?.result.flatMap(result => result.hits)
          .filter(hit => hit.id === this.visibleHitId)
          .flatMap(hit => hit.hsps.map(hsp => ({
            ...hsp,
            hit_id: hit.id,
            definition: hit.definition,
            accession: hit.accession,
            length: hit.length
          }))) || [];
      }
      return this.job?.result.flatMap(result => result.hits).flatMap(hit =>
        hit.hsps.map(hsp => ({
          ...hsp,
          hit_id: hit.id,
          definition: hit.definition,
          accession: hit.accession,
          length: hit.length
        }))
      ) || [];
    }
  },
  async mounted() {
    this.experimentId = this.$route.params.experimentId as string;

    this.job = await getBlastJobApi(this.jobId as string);
    if (this.job) {
      this.editableJobName = this.job.job_name || '';
      this.jobInputs = {
        descriptions: this.job.descriptions,
        alignments: this.job.alignments,
        hitlist_size: this.job.hitlist_size,
        expect: this.job.expect
      };
    }

    this.protein = await getProteinContent(this.job?.protein_id);

    this.jobStatus = await getJobStatus(this.jobId as string);

    this.$nextTick(() => {
      this.visualizeResults();
    });
  },
  methods: {
    async updateJobName() {
      if (this.job && this.editableJobName !== this.job.job_name) {
        try {
          await changeJobName(this.job.job_id, this.editableJobName);
          this.$q.notify({
            type: 'positive',
            message: 'Job name updated successfully.',
          });
        } catch (error) {
          this.editableJobName = this.job.job_name;
          this.$q.notify({
            type: 'negative',
            message: 'Failed to update job name.',
          });
        }
      }
    },
    async updateJobInputs() {
      const setupJobRequest = {
        experiment_id: this.experimentId as string,
        protein_id: this.job?.protein_id as string,
        descriptions: this.jobInputs.descriptions,
        alignments: this.jobInputs.alignments,
        hitlist_size: this.jobInputs.hitlist_size,
        expect: this.jobInputs.expect,
        job_id: this.job?.job_id,
        job_name: this.job?.job_name
      };
      try {
        await setupBlastJob(setupJobRequest as nolabs__application__use_cases__blast__api_models__SetupJobRequest);
        this.$q.notify({
          type: 'positive',
          message: 'Job inputs updated successfully.',
        });
      } catch (error) {
        this.$q.notify({
          type: 'negative',
          message: 'Failed to update job inputs.',
        });
      }
    },
    visualizeResults() {
      if (this.job && this.job.result) {
        const allHSPs = [];
        this.job.result.forEach(result => {
          result.hits.forEach(hit => {
            hit.hsps.forEach(hsp => {
              allHSPs.push({ hitId: hit.id, hitDefinition: hit.definition, accession: hit.accession, ...hsp });
            });
          });
        });
        this.createVisualization('#hit-visualization', allHSPs);
      }
    },
    createVisualization(selector, hsps) {
      const margin = { top: 25, right: 120, bottom: 30, left: 40 }; // Increased right margin for legend
      const width = window.innerWidth - margin.left - margin.right;
      const height = 450 - margin.top - margin.bottom;

      d3.select(selector).selectAll("*").remove(); // Clear the existing graph

      const svg = d3.select(selector)
        .append("svg")
        .attr("width", width + margin.left + margin.right)
        .attr("height", height + margin.top + margin.bottom)
        .append("g")
        .attr("transform", `translate(${margin.left},${margin.top})`);

      const x = d3.scaleLinear()
        .domain([0, d3.max(hsps, d => d.query_to)])
        .range([0, width]);

      const y = d3.scaleBand()
        .domain(hsps.map((hsp, index) => `HSP ${index + 1}`))
        .range([0, height])
        .padding(0.1);

      const color = d3.scaleSequential(d3.interpolateBlues)
        .domain([d3.max(hsps, d => d.score), d3.min(hsps, d => d.score)]); // Light blue for high score to dark blue for low score

      svg.append("g")
        .attr("transform", `translate(0,${height})`)
        .call(d3.axisBottom(x));

      svg.append("g")
        .call(d3.axisLeft(y));

      const bars = svg.selectAll("rect")
        .data(hsps)
        .enter()
        .append("rect")
        .attr("x", d => x(d.query_from))
        .attr("y", (d, i) => y(`HSP ${i + 1}`))
        .attr("width", d => x(d.query_to) - x(d.query_from))
        .attr("height", y.bandwidth())
        .attr("fill", d => color(d.score))
        .attr("stroke", d => d.hitId === this.visibleHitId ? 'green' : 'none')
        .attr("stroke-width", d => d.hitId === this.visibleHitId ? 2 : 0)
        .on("click", (event, d) => this.toggleVisibility(d.hitId))
        .append("title")
        .text(d => `Length: ${d.query_to - d.query_from}, E-value: ${d.evalue}, Score: ${d.score}`);

      svg.selectAll(".bar-text")
        .data(hsps)
        .enter()
        .append("text")
        .attr("class", "bar-text")
        .attr("x", d => x(d.query_from) + 5)
        .attr("y", (d, i) => y(`HSP ${i + 1}`) + y.bandwidth() / 2 + 5)
        .text(d => `${d.accession}`)
        .style("fill", d => d.hitId === this.visibleHitId ? "green" : "#000")
        .style("font-size", "10px");

      // Add a legend for the color scale
      const legendHeight = height;
      const legendWidth = 20;
      const legend = svg.append("g")
        .attr("transform", `translate(${width + 40}, 0)`);

      const legendScale = d3.scaleLinear()
        .domain(color.domain())
        .range([legendHeight, 0]);

      const legendAxis = d3.axisRight(legendScale)
        .ticks(6);

      legend.selectAll("rect")
        .data(d3.range(legendHeight), d => d)
        .enter().append("rect")
        .attr("y", d => d)
        .attr("x", 0)
        .attr("height", 1)
        .attr("width", legendWidth)
        .attr("fill", d => color(legendScale.invert(d)));

      legend.append("g")
        .attr("transform", `translate(${legendWidth}, 0)`)
        .call(legendAxis);

      legend.append("text")
        .attr("x", legendWidth / 2)
        .attr("y", -10)
        .attr("text-anchor", "middle")
        .style("fill", "#fff")
        .style("font-size", "15px")
        .text("Score:");
    },
    toggleVisibility(hitId) {
      this.visibleHitId = this.visibleHitId === hitId ? null : hitId;
      this.visualizeResults(); // Re-render the graph to highlight the selected bar
    },
    async downloadCSV() {
      const hits = this.job?.result.flatMap(result => result.hits).flatMap(hit => hit.hsps) || [];
      const csvContent = "data:text/csv;charset=utf-8,"
        + hits.map(hit => `${hit.num},${hit.bit_score},${hit.score},${hit.evalue},${hit.query_from},${hit.query_to},${hit.hit_from},${hit.hit_to},${hit.query_frame},${hit.hit_frame},${hit.identity},${hit.positive},${hit.gaps},${hit.align_len},${hit.qseq},${hit.hseq},${hit.midline}`).join("\n");

      const encodedUri = encodeURI(csvContent);
      const link = document.createElement("a");
      link.setAttribute("href", encodedUri);
      link.setAttribute("download", "blast_results.csv");
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    },
    copyToClipboard(text) {
      navigator.clipboard.writeText(text).then(() => {
        this.$q.notify({
          type: 'positive',
          message: 'Copied to clipboard!'
        });
      });
    }
  },
  components: {
    QBtn, // Ensure QBtn component is registered
    QTable,
    QTd
  },
});
</script>

<style scoped>
.fasta-content-container {
  width: 100%;
  overflow: hidden;
}

.fasta-content {
  white-space: pre-wrap;
  word-wrap: break-word;
  word-break: break-all;
}

.sequence {
  font-family: monospace;
  background-color: #000000;
  padding: 2px 4px;
  border-radius: 4px;
}

.hit-visualization {
  margin-top: 20px;
}

.full-width {
  width: 100%;
}

.scrollable-cell-content {
  max-width: 250px;
  overflow-x: auto;
  white-space: nowrap;
  display: inline-block;
  vertical-align: middle;
}

.copy-btn {
  margin-left: 8px;
}
</style>
