<script lang="ts">
import {defineComponent} from 'vue'
import {PdbColorSchemas, PdbViews} from "src/components/types";


export default defineComponent({
  name: "PdbViewer",
  computed: {
    PdbViews(): any {
      return PdbViews
    },
    PdbFileName(): string {
      if(this.pdbFileNamePrefix){
        return `${this.pdbFileNamePrefix} - ${this.pdbFile!.name}`;
      }

      return this.pdbFile!.name;
    }
  },
  props: {
    pdbFile: File,
    simulation: Boolean,
    pdbFileNamePrefix: {
      type: String,
      required: false
    }
  },
  data() {
    return {
      selectedRepresentation: PdbViews.default,
      blackBackground: true,
      stage: null
    }
  },
  watch: {
    async selectedRepresentation(newValue, _) {
      await this.render(newValue.key);
    },
    async blackBackground(newValue, _) {
      await this.render(this.selectedRepresentation.key);
    }
  },
  methods: {
    savePdb() {
      if (window.navigator.msSaveOrOpenBlob) {
        window.navigator.msSaveBlob(this.pdbFile!, this.pdbFile!.name);
      } else {
        const elem = window.document.createElement('a');
        elem.href = window.URL.createObjectURL(this.pdbFile!);
        elem.download = this.pdbFile!.name;
        document.body.appendChild(elem);
        elem.click();
        document.body.removeChild(elem);
      }
    },
    async renderSimulation(selectedRepresentation: string) {
      setTimeout(async () => {
        this.$refs.viewport!.innerHTML = '';
        const stage = new NGL.Stage(this.$refs.viewport);
        stage.setParameters({backgroundColor: this.blackBackground ? 'black' : 'white'});
        const pdbComponent = await this.loadFileIntoStage(stage, selectedRepresentation, true);
        this.stage = stage;

        var traj = pdbComponent.addTrajectory().trajectory
        var player = new NGL.TrajectoryPlayer(traj, {
          step: 1,
          timeout: 2000,
          interpolateStep: 50,
          start: 0,
          end: traj.numframes,
          interpolateType: "linear",
          mode: "loop",
          direction: "bounce"
        });
        player.play();
        pdbComponent.autoView();

      }, 100);
    },
    async loadFileIntoStage(stage: any, selectedRepresentation: string, asTrajectory: boolean = false) {
      let pdbComponent: any;
      if (selectedRepresentation === PdbViews.default.key) {
        pdbComponent = await stage.loadFile(this.pdbFile, {defaultRepresentation: true, asTrajectory});
      } else {
        pdbComponent = await stage.loadFile(this.pdbFile, {asTrajectory});
      }

      if (selectedRepresentation !== PdbViews.default.key) {
        pdbComponent.addRepresentation(selectedRepresentation);
      }

      return pdbComponent;
    },
    async renderStatic(selectedRepresentation: string) {
      setTimeout(async () => {
        this.$refs.viewport!.innerHTML = '';
        const stage = new NGL.Stage(this.$refs.viewport);
        stage.setParameters({backgroundColor: this.blackBackground ? 'black' : 'white'});
        this.stage = stage;

        const component = await this.loadFileIntoStage(stage, selectedRepresentation);
        component.autoView();
      }, 50);
    },
    async render(selectedRepresentation: string) {
      if(this.stage){
        this.stage.dispose();
      }

      if (this.simulation) {
        await this.renderSimulation(selectedRepresentation);
      } else {
        await this.renderStatic(selectedRepresentation);
      }
    }
  },
  async mounted() {
    await this.render(this.selectedRepresentation.key);
    new ResizeObserver(async () => {
      setTimeout(async () => {
        await this.render(this.selectedRepresentation.key);
      }, 105);
    }).observe(this.$refs.card.$el);
  }
})
</script>

<template>
  <q-card flat bordered class="my-card">
    <q-card-section ref="card">
      <div class="row items-center">
        <div class="col">
          <div class="text-h6">{{ PdbFileName }}</div>
        </div>
        <div class="col">
          <q-btn color="positive" size="md" outline label="Save pdb" @click="savePdb"/>
        </div>
        <div class="col">
          <q-checkbox v-model="blackBackground" label="Black background" />
        </div>
        <div class="col">
          <q-select v-model="selectedRepresentation"
                    class="q-ml-xs"
                    label-color="positive"
                    color="positive"
                    option-value="key"
                    option-label="title"
                    option-disable="inactive" :options="Object.values(PdbViews)"
                    label="Representation"/>
        </div>
      </div>
    </q-card-section>
    <q-card-section>
      <div ref="viewport" style="width: 100%; min-height: 500px;"></div>
    </q-card-section>
  </q-card>
</template>