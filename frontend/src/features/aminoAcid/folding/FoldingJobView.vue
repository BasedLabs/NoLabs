<template>
  <div v-if="jobLoaded">
    <q-separator></q-separator>
    <q-layout container style="height: 100vh">
      <JobHeader :job-name="job?.name" :on-job-name-change-submit="onJobNameChange">
        <q-btn color="info" size="md" outline label="Folding parameters"
               @click="showInferenceForm = !showInferenceForm"/>
      </JobHeader>
      <q-page-container>
        <div class="row" v-if="jobHasGeneratedData">
          <div class="col-6">
            <div class="q-ma-sm">
              <AminoAcidTable :on-amino-acid-open="setActiveAminoAcid" :rows="job?.aminoAcids"/>
            </div>
          </div>
          <div class="col-6">
            <div class="q-pl-sm q-ma-sm">
              <PdbViewer :pdb-file="activeAminoAcid?.pdbFile" :key="activeAminoAcid?.name"/>
            </div>
          </div>
        </div>
      </q-page-container>
    </q-layout>
    <q-dialog v-model="showInferenceForm" position="right" :maximized="true">
      <q-card>
        <q-card-section class="row items-center q-pb-none">
          <div class="text-h6">Folding parameters</div>
          <q-space/>
          <q-btn icon="close" flat round dense v-close-popup/>
        </q-card-section>
        <q-card-section>
          <AminoAcidInferenceForm :on-submit="onSubmit" :properties="job?.properties"/>
        </q-card-section>
      </q-card>
    </q-dialog>
  </div>
</template>

<script lang="ts">
import {defineComponent} from 'vue'
import {QSpinnerOrbit} from 'quasar';
import JobHeader from "src/components/JobHeader.vue";
import {AminoAcid} from "src/features/aminoAcid/folding/types";
import AminoAcidInferenceForm from "src/features/aminoAcid/AminoAcidInferenceForm.vue";
import AminoAcidTable from "src/features/aminoAcid/AminoAcidTable.vue";
import useFoldingStore from "src/features/aminoAcid/folding/storage";
import PdbViewer from "src/components/PdbViewer.vue";
import {Job} from "src/features/aminoAcid/types";


export default defineComponent({
  name: 'LocalisationJobView',
  data() {
    const store = useFoldingStore();

    return {
      job: null as Job<AminoAcid>,
      showInferenceForm: false,
      store,
      activeAminoAcid: null as AminoAcid | null | undefined
    }
  },
  computed: {
    jobLoaded(): boolean {
      return this.job !== null;
    },
    jobHasGeneratedData(): boolean {
      return this.jobLoaded && this.job!.aminoAcids.length > 0;
    },
    aminoAcidRows() {
      return this.job?.aminoAcids;
    },
  },
  methods: {
    setActiveAminoAcid(aminoAcidName: string): void {
      this.activeAminoAcid = this.job?.aminoAcids.find(x => x.name === aminoAcidName);
    },
    setJob(job: Job<AminoAcid> | null) {
      if (job !== null) {
        this.job = job;

        if (job.aminoAcids.length > 0) {
          this.setActiveAminoAcid(job.aminoAcids[0]!.name);
        }
      }
    },
    async onJobNameChange(newJobName: string) {
      await this.store.changeJobName(this.job?.id as string, newJobName);
      this.job!.name = newJobName;
    },
    async onSubmit(data: { fastas: Array<File> }) {
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: 'Running AI models. This can take a couple of minutes'
      });

      const response = await this.store.inference({
        jobId: this.job!.id,
        jobName: this.job!.name,
        fastas: data.fastas
      });

      this.setJob(response.job);

      this.showInferenceForm = false;

      this.$q.loading.hide();
    },
  },
  async mounted() {
    const jobId = this.$route.params.jobId as string;

    this.$q.loading.show({
      spinner: QSpinnerOrbit,
      message: `Job ${jobId}`
    });

    const response = await this.store.getJob(jobId);

    this.setJob(response.job);

    this.$q.loading.hide();

    if (!this.jobHasGeneratedData) {
      this.showInferenceForm = true;
    }
  },
  components: {
    PdbViewer,
    AminoAcidTable,
    JobHeader,
    AminoAcidInferenceForm
  }
})
</script>
