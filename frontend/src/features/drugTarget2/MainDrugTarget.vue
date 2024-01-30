<template>
  <div class="q-pa-md">
    <q-stepper
        v-model="step"
        ref="stepper"
        color="positive"
        animated
    >
      <q-step
          :name="1"
          title="Upload targets"
          icon="cloud-upload"
          color="positive"
          :done="step > 1"
      >
        Upload target files in .fasta format
        <q-file
            v-model="files"
            accept=".fasta"
            label="Pick fasta files"
            filled
            multiple
            style="max-width: 300px"
        />
      </q-step>

      <q-step
          :name="2"
          title="Add ligands"
          color="positive"
          icon="cloud"
          :done="step > 2"
      >
        <q-table
            flat bordered
            title="Target-ligand pairs"
            :rows="targetLigandRows"
            :columns="targetLigandColumns"
            row-key="name"
        >
          <template v-slot:top-right>
            <q-btn size="lg" color="positive"
                   :disable="targetLigandRows.map(x => x.ligands.length > 0).filter(x => x) == 0" style="width: 100%"
                   class="q-mb-xs" label="Run all jobs"/>
          </template>
          <template v-slot:header="props">
            <q-tr :props="props">
              <q-th auto-width/>
              <q-th auto-width/>
              <q-th
                  v-for="col in props.cols"
                  :key="col.name"
                  :props="props"
              >
                {{ col.label }}
              </q-th>
              <q-th auto-width/>
            </q-tr>
          </template>

          <template v-slot:body="props">
            <q-tr :props="props">
              <q-td auto-width>
                <q-btn size="md" flat color="positive" :disable="props.row.ligands.length == 0" style="width: 100%"
                       class="q-mb-xs" label="Run single job"/>
              </q-td>
              <q-td auto-width>
                <q-file
                    v-model="props.row.ligands"
                    accept=".sdf"
                    label="Pick ligand files"
                    filled
                    multiple
                    style="width: 200px"
                />
              </q-td>
              <q-td
                  v-for="col in props.cols"
                  :key="col.name"
                  :props="props"
              >
                {{ col.value }}
              </q-td>
            </q-tr>
            <q-tr v-show="props.expand" :props="props">
              <q-td colspan="100%">
                <div class="text-left">This is expand slot for row above: {{ props.row.name }}.</div>
              </q-td>
            </q-tr>
          </template>

        </q-table>
      </q-step>

      <template v-slot:navigation>
        <q-stepper-navigation>
          <q-btn @click="$refs.stepper.next()" color="positive" :label="step === 4 ? 'Finish' : 'Continue'"/>
          <q-btn v-if="step > 1" flat color="positive" @click="$refs.stepper.previous()" label="Back" class="q-ml-sm"/>
        </q-stepper-navigation>
      </template>
    </q-stepper>
  </div>
</template>

<script lang="ts">
import {defineComponent} from "vue";

export default defineComponent({
  name: "MainDrugTarget",
  data() {
    return {
      step: 1,
      files: [] as Array<File>,
      targetLigandRows: [] as Array<{ targetName: string, ligands: Array<File> }>,
      targetLigandColumns: [
        {
          name: 'targetName',
          required: true,
          label: 'Target',
          align: 'left',
          field: (row: { targetName: string, ligands: Array<File> }) => row.targetName
        },
        {
          name: 'ligandsNames',
          required: true,
          label: 'Ligands',
          align: 'left',
          field: (row: { targetName: string, ligands: Array<File> }) => row.ligands.map(x => x.name).join(', ')
        }
      ]
    };
  },
  watch: {
    files(newValue, _) {
      let result = [];
      for (let i = 0; i < newValue.length; i++) {
        const file = newValue[i];
        result.push({
          targetName: file.name,
          ligands: []
        })
      }
      this.targetLigandRows = result;
    }
  }
})
</script>
