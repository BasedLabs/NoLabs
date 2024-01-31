<template>
  <q-card v-if="isDocking">
    <q-card-section>
      <q-stepper v-model="dockingStep" vertical color="info" animated>
        <q-step
          color="info"
          :name="1"
          title="Predict MSA"
          :done="dockingStep > 1"
        >
          <q-linear-progress :value="progress[0]" striped color="blue" />
        </q-step>
        <q-step
          color="info"
          :name="2"
          title="Get Folded Structure"
          :done="dockingStep > 2"
        >
          <q-linear-progress :value="progress[1]" striped color="green" />
        </q-step>
        <q-step
          color="info"
          :name="3"
          title="Get Binding Pocket"
          :done="dockingStep > 3"
        >
          <q-linear-progress :value="progress[2]" striped color="orange" />
        </q-step>
        <q-step
          :name="4"
          color="info"
          title="Run Docking"
          :done="dockingStep > 4"
        >
          <q-linear-progress :value="progress[3]" striped color="red" />
          <q-btn label="View Results" @click="viewResults" />
          <q-card v-if="showResults">
            <q-card-section>
              <div ref="resultsViewContainer" class="results-viewer"></div>
            </q-card-section>
            <q-card-actions align="right">
              <q-btn
                flat
                label="Download Results"
                color="primary"
                @click="downloadResults"
              />
              <q-btn flat label="Close" color="primary" v-close-popup />
            </q-card-actions>
          </q-card>
        </q-step>
      </q-stepper>
    </q-card-section>
  </q-card>
</template>

<script>
import * as NGL from "ngl";
export default {
  name: "DockingDetail",
  components: {},
  data() {
    return {
      loading: true,
      isDocking: false,
      dockingStep: 0,
      progress: [0, 0, 0, 0],
      dockingTimers: [10, 1, 1, 30], // Seconds for each step
      showResults: false,
      resultsViewer: null,
    };
  },
  methods: {
    startProgress(stepIndex) {
      const duration = this.dockingTimers[stepIndex] * 1000; // Convert seconds to milliseconds
      const interval = 100; // Update interval in milliseconds
      let elapsed = 0;

      const timer = setInterval(() => {
        elapsed += interval;
        this.progress[stepIndex] = elapsed / duration;

        if (elapsed >= duration) {
          clearInterval(timer);
          if (stepIndex < this.dockingTimers.length - 1) {
            this.dockingStep += 1; // Move to next step
            this.startProgress(stepIndex + 1); // Start next step's progress
          }
        }
      }, interval);
    },
    viewResults() {
      this.showResults = true;
      this.$nextTick(() => {
        if (!this.resultsViewer) {
          this.resultsViewer = new NGL.Stage(this.$refs.resultsViewContainer, {
            backgroundColor: "black",
          });

          this.resultsViewer
            .loadFile("/7NB4_pred_ligand.sdf")
            .then((component) => {
              component.addRepresentation("ball+stick");
              component.autoView();
            });

          this.resultsViewer
            .loadFile("/7NB4_pred_protein.pdb")
            .then((component) => {
              component.addRepresentation("cartoon", { color: "sstruc" });
              component.autoView();
            });
        }
      });
    },
  },
  mounted() {
    this.isDocking = true;
    this.dockingStep = 1;
    this.startProgress(0); // Start progress for the first step
  },
};
</script>
