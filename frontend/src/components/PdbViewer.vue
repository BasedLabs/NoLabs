<script lang="ts">
import { defineComponent } from "vue";
import { PdbViews } from "src/components/types";

export default defineComponent({
  name: "PdbViewer",
  computed: {
    PdbViews(): any {
      return PdbViews;
    },
  },
  props: {
    pdbFile: File,
  },
  data() {
    return {
      selectedRepresentation: PdbViews.default,
    };
  },
  watch: {
    async selectedRepresentation(newValue, _) {
      await this.render(newValue.key);
    },
  },
  methods: {
    savePdb() {
      if (window.navigator.msSaveOrOpenBlob) {
        window.navigator.msSaveBlob(this.pdbFile!, this.pdbFile!.name);
      } else {
        const elem = window.document.createElement("a");
        elem.href = window.URL.createObjectURL(this.pdbFile!);
        elem.download = this.pdbFile!.name;
        document.body.appendChild(elem);
        elem.click();
        document.body.removeChild(elem);
      }
    },
    async render(selectedRepresentation: string) {
      setTimeout(async () => {
        this.$refs.viewport.innerHTML = "";
        const stage = new NGL.Stage(this.$refs.viewport);
        stage.setParameters({ backgroundColor: "black" });
        if (selectedRepresentation === PdbViews.default.key) {
          await stage.loadFile(this.pdbFile, { defaultRepresentation: true });
        } else {
          const pdbComponent = await stage.loadFile(this.pdbFile);
          pdbComponent.addRepresentation(selectedRepresentation);
        }
      }, 100);
    },
  },
  async mounted() {
    await this.render(this.selectedRepresentation.key);
  },
});
</script>

<template>
  <q-card flat bordered class="my-card">
    <q-card-section>
      <div class="row items-center">
        <div class="col">
          <div class="text-h6" align="-">{{ pdbFile.name }}</div>
        </div>
        <div class="col">
          <q-btn
            color="positive"
            size="md"
            outline
            label="Save pdb"
            @click="savePdb"
          />
        </div>
        <div class="col">
          <q-select
            v-model="selectedRepresentation"
            label-color="positive"
            color="positive"
            option-value="key"
            option-label="title"
            option-disable="inactive"
            :options="Object.values(PdbViews)"
            label="Representation"
          />
        </div>
      </div>
    </q-card-section>
    <q-card-section>
      <div ref="viewport" style="width: 100%; min-height: 400px"></div>
    </q-card-section>
  </q-card>
</template>
