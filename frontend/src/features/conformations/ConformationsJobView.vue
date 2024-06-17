<template>
  <div v-if="jobLoaded">
    <q-separator></q-separator>
    <q-layout container style="height: 100vh">
      <JobHeader :job-name="job!.name!" :on-job-name-change-submit="onJobNameChange">
        <q-btn color="info" size="md" outline label="Simulation parameters"
               @click="showInferenceForm = !showInferenceForm"/>
      </JobHeader>
      <q-page-container>
        <div class="row">
          <div :class="tiles.one.current"
               style="transition: all .1s linear;">
            <div class="q-ma-sm">
              <q-btn size="xs" flat color="info" style="width: 100%" class="q-mb-xs" label="EXPAND"
                     @click="expandTile('one');"/>
              <PdbViewer v-if="jobHasInputData" :pdb-file="job!.properties!.pdbFile"
                         :key="job!.properties.pdbFile!.name"/>
            </div>
          </div>
          <div :class="tiles.two.current"
               style="transition: all .1s linear;">
            <div class="q-pl-sm q-ma-sm">
              <q-btn size="xs" flat color="info" style="width: 100%" class="q-mb-xs" label="EXPAND"
                     @click="expandTile('two');"/>
              <TimelineView :timeline="job!.timeline!"/>
            </div>
          </div>
          <div :class="tiles.three.current" style="transition: all .1s linear;">
            <div class="q-mt-sm q-mb-sm q-mr-sm">
              <q-btn flat size="xs" color="info" style="width: 100%" class="q-mb-xs" label="EXPAND"
                     @click="expandTile('three');"/>
              <h6 v-if="simulationCriticalError" class="q-ma-md">
                <q-icon name="warning" color="warning" size="4rem"/>
                Simulations error. Fix errors (check them in table) and restart simulation.
              </h6>
              <PdbViewer v-if="!simulationCriticalError && jobHasGeneratedData"
                         :key="job?.pdbContent?.name"
                         :pdb-file="job?.pdbContent" :simulation="true"
                         file-name-prefix="Simulation"/>
            </div>
          </div>
        </div>
      </q-page-container>
    </q-layout>
    <q-dialog v-model="showInferenceForm" position="right" :maximized="true">
      <q-card>
        <q-card-section class="row items-center q-pb-none">
          <div class="text-h6">Simulation parameters</div>
          <q-space/>
          <q-btn icon="close" flat round dense v-close-popup/>
        </q-card-section>
        <q-card-section>
          <InferenceFormView :on-submit="onSubmit" :properties="job!.properties!" :save-parameters="saveParameters"/>
        </q-card-section>
      </q-card>
    </q-dialog>
  </div>
</template>

<script lang="ts">
import {defineComponent, PropType} from 'vue'
import {QVueGlobals, QSpinnerOrbit} from 'quasar';
import PdbViewer from "src/components/PdbViewer.vue";
import InferenceFormView from "src/features/conformations/InferenceFormView.vue";
import useConformationsStore from "src/features/conformations/storage";
import {Job, JobProperties} from "src/features/conformations/types";
import TimelineView from "src/components/Timeline.vue";
import JobHeader from "src/components/JobHeader.vue";


export default defineComponent({
  name: 'ConformationJobView',
  data() {
    const store = useConformationsStore();

    return {
      job: null as Job,
      showInferenceForm: false,
      store,
      tiles: {
        one: {
          current: 'col-4',
          hover: 'col-6',
          leave: 'col-4',
          otherHover: 'col-3'
        },
        two: {
          current: 'col-4',
          hover: 'col-6',
          leave: 'col-4',
          otherHover: 'col-3'
        },
        three: {
          current: 'col-4',
          hover: 'col-6',
          leave: 'col-4',
          otherHover: 'col-3'
        }
      } as { [index: string]: { current: String, hover: String, leave: String, otherHover: String } }
    }
  },
  computed: {
    jobLoaded(): boolean {
      return this.job !== null;
    },
    jobHasInputData(): boolean {
      return this.jobLoaded && this.job!.properties.pdbFile != null;
    },
    jobHasGeneratedData(): boolean {
      return this.jobLoaded && this.job!.pdbContent != null && this.job!.pdbContent!.size > 4;
    },
    simulationCriticalError(): boolean {
      return this.jobLoaded && this.job!.timeline.length > 0 && (this.job?.pdbContent == null || this.job?.pdbContent.size <= 4);
    }
  },
  methods: {
    expandTile(index: string) {
      if (this.tiles[index].current === this.tiles[index].hover) {
        for (const key in this.tiles) {
          this.tiles[key].current = this.tiles[key].leave;
        }
        return;
      }

      for (const key in this.tiles) {
        if (key == index) {
          this.tiles[key].current = this.tiles[key].hover;
        } else {
          this.tiles[key].current = this.tiles[key].otherHover;
        }
      }
    },
    async saveParameters(properties: JobProperties){
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: 'Saving parameters'
      });

      await this.store.saveParameters({
        jobId: this.job!.id,
        jobName: this.job!.name,
        pdbFile: properties.pdbFile!,
        totalFrames: properties.totalFrames!,
        temperatureK: properties.temperatureK!,
        takeFrameEvery: properties.takeFrameEvery!,
        stepSize: properties.stepSize!,
        replaceNonStandardResidues: properties.replaceNonStandardResidues!,
        addMissingAtoms: properties.addMissingAtoms!,
        addMissingHydrogens: properties.addMissingHydrogens!,
        frictionCoeff: properties.frictionCoeff!,
        ignoreMissingAtoms: properties.ignoreMissingAtoms!,
        integrator: properties.integrator!,
      });

      this.showInferenceForm = false;

      this.$q.loading.hide();
    },
    async onSubmit(properties: JobProperties) {
      this.$q.loading.show({
        spinner: QSpinnerOrbit,
        message: 'Running computations'
      });

      const response = await this.store.inference({
        jobId: this.job!.id,
        jobName: this.job!.name,
        pdbFile: properties.pdbFile!,
        totalFrames: properties.totalFrames!,
        temperatureK: properties.temperatureK!,
        takeFrameEvery: properties.takeFrameEvery!,
        stepSize: properties.stepSize!,
        replaceNonStandardResidues: properties.replaceNonStandardResidues!,
        addMissingAtoms: properties.addMissingAtoms!,
        addMissingHydrogens: properties.addMissingHydrogens!,
        frictionCoeff: properties.frictionCoeff!,
        ignoreMissingAtoms: properties.ignoreMissingAtoms!,
        integrator: properties.integrator!,
      });

      if (response.job !== null) {
        this.job = response.job;
      }

      this.showInferenceForm = false;

      this.$q.loading.hide();
    },
    async onJobNameChange(newJobName: string) {
      await this.store.changeJobName(this.job?.id as string, newJobName);
      this.job!.name = newJobName;
    }
  },
  async mounted() {
    const jobId = this.$route.params.jobId as string;

    this.$q.loading.show({
      spinner: QSpinnerOrbit,
      message: `Job ${jobId}`
    });

    const response = await this.store.getJob(jobId);

    if (response.job !== null) {
      this.job = response.job;
    }

    this.$q.loading.hide();

    if (!this.jobHasInputData) {
      this.showInferenceForm = true;
    }
  },
  components: {
    JobHeader,
    TimelineView,
    PdbViewer,
    InferenceFormView
  }
})
</script>
