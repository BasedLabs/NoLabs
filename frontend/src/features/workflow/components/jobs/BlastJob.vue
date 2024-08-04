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
                  <q-btn v-if="!jobStatus.running" @click="startJob" color="info" label="Start" />
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
              <q-item-section>{{ job?.descriptions }}</q-item-section>
            </q-item>
            <q-item>
              <q-item-section>
                <q-item-label>Alignments</q-item-label>
              </q-item-section>
              <q-item-section>{{ job?.alignments }}</q-item-section>
            </q-item>
            <q-item>
              <q-item-section>
                <q-item-label>Hitlist Size</q-item-label>
              </q-item-section>
              <q-item-section>{{ job?.hitlist_size }}</q-item-section>
            </q-item>
            <q-item>
              <q-item-section>
                <q-item-label>Expect</q-item-label>
              </q-item-section>
              <q-item-section>{{ job?.expect }}</q-item-section>
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
          <div v-if="jobHasGeneratedData">
            <div id="hit-visualization" class="hit-visualization"></div>
            <div v-for="result in job?.result" :key="result.protein_id">
              <div v-for="hit in result.hits" :key="hit.id" class="hit-section" v-show="isVisible(hit.id)">
                <h6>Hit ID: {{ hit.id }}</h6>
                <p>Definition: {{ hit.definition }}</p>
                <p>Accession: {{ hit.accession }}</p>
                <p>Length: {{ hit.length }}</p>
                <div v-for="hsp in hit.hsps" :key="hsp.num" class="hsp-section">
                  <p>HSP Number: {{ hsp.num }}</p>
                  <p>Bit Score: {{ hsp.bit_score }}</p>
                  <p>Score: {{ hsp.score }}</p>
                  <p>E-value: {{ hsp.evalue }}</p>
                  <p>Query From: {{ hsp.query_from }}</p>
                  <p>Query To: {{ hsp.query_to }}</p>
                  <p>Hit From: {{ hsp.hit_from }}</p>
                  <p>Hit To: {{ hsp.hit_to }}</p>
                  <p>Query Frame: {{ hsp.query_frame }}</p>
                  <p>Hit Frame: {{ hsp.hit_frame }}</p>
                  <p>Identity: {{ hsp.identity }}</p>
                  <p>Positive: {{ hsp.positive }}</p>
                  <p>Gaps: {{ hsp.gaps }}</p>
                  <p>Alignment Length: {{ hsp.align_len }}</p>
                  <p>Query Sequence: <span class="sequence">{{ hsp.qseq }}</span></p>
                  <p>Hit Sequence: <span class="sequence">{{ hsp.hseq }}</span></p>
                  <p>Midline: <span class="sequence">{{ hsp.midline }}</span></p>
                </div>
              </div>
            </div>
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
import { QSpinner, QInput, QBtn } from 'quasar';
import * as d3 from 'd3';
import {
  nolabs__application__use_cases__blast__api_models__JobResponse,
  ProteinContentResponse,
  nolabs__application__use_cases__blast__api_models__GetJobStatusResponse,
  nolabs__application__use_cases__blast__api_models__SetupJobRequest
} from "src/refinedApi/client";
import {
  getBlastJobApi,
  getProteinContent,
  getBlastJobStatus,
  changeJobName,
  setupBlastJob,
  startBlastJob
} from "src/features/workflow/refinedApi";

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
  },
  async mounted() {
    this.experimentId = this.$route.params.experimentId as string;

    this.job = await getBlastJobApi(this.jobId as string);
    if (this.job) {
      this.editableJobName = this.job.job_name || '';
    }

    this.protein = await getProteinContent(this.job?.protein_id);

    this.jobStatus = await getBlastJobStatus(this.jobId as string);

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
    async startJob() {
      try {
        await startBlastJob(this.jobId as string);
        this.$q.notify({
          type: 'positive',
          message: 'Job started successfully.',
        });
        this.jobStatus = await getBlastJobStatus(this.jobId as string);
      } catch (error) {
        this.$q.notify({
          type: 'negative',
          message: 'Failed to start the job.',
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
      const margin = { top: 10, right: 120, bottom: 30, left: 40 }; // Increased right margin for legend
      const width = 800 - margin.left - margin.right;
      const height = 400 - margin.top - margin.bottom;

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

      const color = d3.scaleSequential(d3.interpolateCool)
        .domain([d3.min(hsps, d => d.score), d3.max(hsps, d => d.score)]);

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
        .on("click", (event, d) => this.toggleVisibility(d.hitId));

      svg.selectAll(".bar-text")
        .data(hsps)
        .enter()
        .append("text")
        .attr("class", "bar-text")
        .attr("x", d => x(d.query_from) + 5)
        .attr("y", (d, i) => y(`HSP ${i + 1}`) + y.bandwidth() / 2 + 5)
        .text(d => `${d.accession}`)
        .style("fill", "#000")
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
        .style("font-size", "12px")
        .text("Score:");
    },
    toggleVisibility(hitId) {
      this.visibleHitId = this.visibleHitId === hitId ? null : hitId;
    },
    isVisible(hitId) {
      return this.visibleHitId === hitId;
    }
  },
  components: {
    QBtn, // Ensure QBtn component is registered
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

.hit-section {
  border: 1px solid #00e1ff;
  padding: 10px;
  margin-bottom: 10px;
}

.hsp-section {
  border: 1px solid #00e1ff;
  padding: 10px;
  margin-bottom: 10px;
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
</style>
