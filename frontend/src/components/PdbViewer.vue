<script lang="ts">
import {defineComponent, PropType} from 'vue'
import {PdbViews} from "src/components/types";
import {ExperimentProperties} from "../features/proteinDesign/types";


interface OnAtomClick {
  (stage: any, coordinates: {x: number, y: number, z: number}): void;
}

interface OnRender {
  (stage): void;
}


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
    onRender: {
      type: Function as PropType<OnRender>,
      required: false
    },
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
    },
    onAtomClick: {
      type: Function as PropType<OnAtomClick>,
      required: false
    }
  },
  data() {
    return {
      selectedRepresentation: PdbViews.default,
      blackBackground: true,
      stage: null,
      lastClickedAtomPosition: null,
      loading: false
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
      if (!this.pdbFile?.size) {
        return;
      }

      setTimeout(async () => {
        this.loading = true;
        if (this.$refs.viewport) {
          this.$refs.viewport!.innerHTML = '';
        }
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
        this.loading = false;

        if(this.onRender){
            this.onRender!(stage);
        }
      }, 500);
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

      this.$refs.viewport.addEventListener('mouseover', this.handleMouseover)
      return component;
    },
    async renderStatic(selectedRepresentation: string) {
      if (!this.pdbFile?.size) {
        return;
      }

      setTimeout(async () => {
        this.loading = true;
        if (this.$refs.viewport) {
          this.$refs.viewport!.innerHTML = '';
        }
        const stage = new NGL.Stage(this.$refs.viewport);
        stage.setParameters({backgroundColor: this.blackBackground ? 'black' : 'white'});
        this.$options.stage = stage;

        stage.signals.clicked.add((pickingProxy) => {
          // Check if the pickingProxy is valid and represents an atom
          if (pickingProxy && pickingProxy.atom) {
            const vector = pickingProxy.atom.positionToVector3();
            this.lastClickedAtomPosition = vector

            if(this.onAtomClick){
              this.onAtomClick(stage, vector);
            }
          }
        });

        const component = await this.loadFileIntoStage(stage, selectedRepresentation);
        component.autoView();
        this.loading = false;

        if(this.onRender){
            this.onRender!(stage);
        }
      }, 500);
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
        <div class="col-1">
          <q-btn color="info" icon="help" outline round size="sm">
            <q-tooltip :offset="[10, 10]">
              Click on atom to show its position.
              Right click on atom 1 + double right click on atom 2 to show distance in Å.
            </q-tooltip>
          </q-btn>
        </div>
        <div class="col-1" v-if="pdbFile != null">
          <q-btn color="info" size="sm" outline icon="file_download" @click="savePdb"/>
        </div>
        <div class="col">
          <div class="text-h6">{{ FileName }}</div>
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
        <div class="col" v-if="lastClickedAtomPosition">
          Atom position
          {{ lastClickedAtomPosition.x.toFixed(2) }},{{
            lastClickedAtomPosition.y.toFixed(2)
          }},{{ lastClickedAtomPosition.z.toFixed(2) }}
        </div>
      </div>
    </q-card-section>
    <q-card-section>
      <div ref="viewport" style="width: 100%; min-height: 500px;"></div>
      <q-inner-loading
        :showing="loading"
        label="Loading"
      />
    </q-card-section>
  </q-card>
</template>
