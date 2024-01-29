<script lang="ts">
import {defineComponent, PropType} from 'vue'
import {AminoAcid} from "src/features/aminoAcid/geneOntology/types";

export default defineComponent({
  name: "AminoAcidTable",
  props: {
    rows: {
      type: Object as PropType<{name: string, sequence: string}>,
      required: true
    },
    onAminoAcidOpen: {
      type: Function as PropType<(aminoAcidName: string) => void>,
      required: true
    }
  },
  methods: {
    aminoAcidRowClick(aminoAcid: {name: string, sequence: string}){
      this.selected = [aminoAcid];
      this.onAminoAcidOpen(aminoAcid.name);
    }
  },
  data() {
    return {
      selected: [] as Array<{name: string, sequence: string}>,
      columns: [
        {
          name: 'name',
          label: 'Name',
          align: 'left',
          sortable: 'false',
          field: (row: AminoAcid) => row.name,
        }, {
          name: 'sequence',
          label: 'Sequence',
          align: 'left',
          sortable: true,
          field: (row: AminoAcid) => row.sequence,
        }
      ]
    }
  }
})
</script>

<template>
  <q-table
      title="Amino acids"
      :rows="rows"
      :columns="columns"
      row-key="name"
      v-model:selected="selected"
      :rows-per-page-options="[3, 5, 10]"
  >
    <template v-slot:body="props">
      <q-tr :props="props" @click="aminoAcidRowClick(props.row)">
        <q-td
            auto-width
            v-for="col in props.cols"
            :key="col.name"
            :props="props"
            class="hover-finger"
        >
          {{ col.value }}
        </q-td>
      </q-tr>
    </template>
  </q-table>
</template>