<script lang="ts">
import {defineComponent, PropType} from 'vue'
import {PdbViews} from "src/components/types";


export default defineComponent({
  name: "PdbViewer",
  stage: null,
  computed: {
    PdbViews(): any {
      return PdbViews
    },
    FileName(): string {
      const name = this.pdbFile != null ? this.pdbFile!.name : this.sdfFile!.name;
      if (this.fileNamePrefix) {
        return `${this.fileNamePrefix} - ${name}`;
      }

      return name;
    }
  },
  props: {
    pdbFile: {
      type: Object as PropType<File | null>,
      required: false
    },
    sdfFile: {
      type: Object as PropType<File | null>,
      required: false
    },
    pocketIds: {
      type: [] as PropType<Array<number>>,
      required: false
    },
    simulation: Boolean,
    fileNamePrefix: {
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
        this.$options.stage = stage;

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
        pdbComponent.addRepresentation("unitcell");
        pdbComponent.autoView();

      }, 100);
    },
    async loadFileIntoStage(stage: any, selectedRepresentation: string, asTrajectory: boolean = false) {
      let component: any;
      if (this.pdbFile != null) {
        if (selectedRepresentation === PdbViews.default.key) {
          component = await stage.loadFile(this.pdbFile, {defaultRepresentation: true, asTrajectory});
        } else {
          component = await stage.loadFile(this.pdbFile, {asTrajectory});
        }
        if (selectedRepresentation !== PdbViews.default.key) {
          component.addRepresentation(selectedRepresentation);
        }
      }
      if (this.sdfFile != null) {
        if (selectedRepresentation === PdbViews.default.key) {
          component = await stage.loadFile(this.sdfFile, {defaultRepresentation: true});
        } else {
          component = await stage.loadFile(this.sdfFile, {asTrajectory});
        }
        if (selectedRepresentation !== PdbViews.default.key) {
          component.addRepresentation(selectedRepresentation);
        }
      }
      if (this.pocketIds && this.pocketIds.length > 0) {
        const selectionString = this.pocketIds
          .map((id) => (id + 1).toString()).join(" or ");
        component.addRepresentation("ball+stick", {
          sele: selectionString,
          color: "blue",
        });
      }
      component.addRepresentation("unitcell");
      this.$refs.viewport.addEventListener('mouseover', this.handleMouseover)
      return component;
    },
    async renderStatic(selectedRepresentation: string) {
      setTimeout(async () => {
        this.$refs.viewport!.innerHTML = '';
        const stage = new NGL.Stage(this.$refs.viewport);
        stage.setParameters({backgroundColor: this.blackBackground ? 'black' : 'white'});
        this.$options.stage = stage;

        const component = await this.loadFileIntoStage(stage, selectedRepresentation);
        component.autoView();
      }, 100);
    },
    async render(selectedRepresentation: string) {
      if (this.$options.stage) {
        this.$options.stage.dispose();
      }

      if (this.simulation) {
        await this.renderSimulation(selectedRepresentation);
      } else {
        await this.renderStatic(selectedRepresentation);
      }
    },
    handleMouseover(event: any) {
      // Check if the stage is available
      if (this.$options.stage) {
        console.log(this.$options.stage);
        // Get mouse position relative to the molecular structure
        const mousePosition = this.$options.stage.viewerControls.getMousePosition(event);
        if (mousePosition) {
          // Extract X, Y, Z coordinates
          const { x, y, z } = mousePosition;
          // Update UI with the coordinates
          this.updateMouseCoordinates({ x, y, z });
        }
      }
    }
  },
  async mounted() {
    if (this.pdbFile == null && this.sdfFile == null) {
      console.error('Specify pdbFile or sdfFile, or both');
      return;
    }

    if (this.sdfFile != null && this.simulation) {
      console.error('Cannot run simulation for sdf file');
      return;
    }

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
          <div class="text-h6">{{ FileName }}</div>
        </div>
        <div class="col" v-if="pdbFile != null">
          <q-btn color="info" size="md" outline label="Save pdb" @click="savePdb"/>
        </div>
        <div class="col">
          <q-checkbox v-model="blackBackground" label="Black background"/>
        </div>
        <div class="col">
          <q-select v-model="selectedRepresentation"
                    class="q-ml-xs"
                    label-color="info"
                    color="info"
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
