<script lang="ts">
import {defineComponent, PropType} from 'vue'
import {Timeline} from "src/components/types";

export default defineComponent({
  name: "TimelineView",
  props: {
    timeline: {
      type: Array<Timeline>,
      required: true
    }
  },
  computed: {
    timelineColumns(): any {
      return [
        {
          name: 'createdAt',
          required: true,
          label: 'Created at',
          align: 'left',
          sortable: true,
          field: (row: Timeline) => row.createdAt,
        },
        {
          name: 'message',
          required: true,
          label: 'Message',
          align: 'left',
          sortable: false,
          field: (row: Timeline) => row.message
        }
      ]
    }
  }
})
</script>

<template>
  <q-table
      flat bordered
      title="Timeline and errors"
      :rows="timeline"
      :columns="timelineColumns"
      :wrap-cells="true"
      row-key="name"
      :pagination="{
        rowsPerPage: 10,
        sortBy: 'createdAt',
        descending: true,
      }"
  >

    <template v-slot:header="props">
      <q-tr :props="props">
        <q-th auto-width/>
        <q-th
            v-for="col in props.cols"
            :key="col.name"
            :props="props"
        >
          {{ col.label }}
        </q-th>
      </q-tr>
    </template>

    <template v-slot:body="props">
      <q-tr :props="props">
        <q-td auto-width>
          <q-btn v-if="props.row.error !== null && props.row.error !== ''" size="sm" color="negative" dense
                 @click="props.row.expand = !props.row.expand" :icon="props.row.expand ? 'remove' : 'add'"/>
        </q-td>
        <q-td
            v-for="col in props.cols"
            :key="col.name"
            :props="props"
        >
          {{ col.value }}
        </q-td>
      </q-tr>
      <q-tr v-show="props.row.expand" :props="props">
        <q-td colspan="100%">
          <div class="text-left">{{ props.row.error }}</div>
        </q-td>
      </q-tr>
    </template>
  </q-table>
</template>