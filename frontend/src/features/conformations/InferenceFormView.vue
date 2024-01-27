<script lang="ts">
import { defineComponent, PropType } from "vue";
import { IntegratorsRequest } from "src/api/client";
import { ExperimentProperties } from "src/features/conformations/types";

interface OnSubmitType {
  (inputs: ExperimentProperties): Promise<void>;
}

export default defineComponent({
  name: "InferenceFormView",
  props: {
    onSubmit: {
      type: Function as PropType<OnSubmitType>,
      required: true,
    },
    properties: {
      type: {} as PropType<ExperimentProperties>,
      required: true,
    },
  },
  computed: {
    IntegratorsValues() {
      return Object.values(IntegratorsRequest);
    },
    pdbFileRule() {
      return this.pdbFile !== null;
    },
  },
  data(): ExperimentProperties {
    return {
      pdbFile: this.properties!.pdbFile,
      totalFrames: this.properties!.totalFrames,
      temperatureK: this.properties!.temperatureK,
      takeFrameEvery: this.properties!.takeFrameEvery,
      stepSize: this.properties!.stepSize,
      replaceNonStandardResidues: this.properties!.replaceNonStandardResidues,
      addMissingAtoms: this.properties!.addMissingAtoms,
      addMissingHydrogens: this.properties!.addMissingHydrogens,
      frictionCoeff: this.properties!.frictionCoeff,
      ignoreMissingAtoms: this.properties!.ignoreMissingAtoms,
      integrator: this.properties!.integrator,
    };
  },
  methods: {
    async _onSubmit() {
      await this.onSubmit!(this.$data);
    },
  },
});
</script>

<template>
  <q-form @submit="_onSubmit" class="q-gutter-md">
    <q-input
      filled
      v-model="totalFrames"
      label="Total frames"
      lazy-rules
      type="number"
      :rules="[(val) => (val && val.length > 0) || 'Please type something']"
    >
    </q-input>
    <q-input
      filled
      v-model="temperatureK"
      label="Temperature"
      type="number"
      :rules="[(val) => (val && val.length > 0) || 'Please type something']"
    >
      <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
        Temperature of the system
      </q-tooltip>
    </q-input>
    <q-input
      filled
      type="number"
      v-model="takeFrameEvery"
      label="Frame snapshot rate"
      lazy-rules
      :rules="[(val) => (val && val > 0) || 'Please type something']"
    >
      <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
        Take frame every N iterations
      </q-tooltip>
    </q-input>
    <q-input
      filled
      type="number"
      v-model="stepSize"
      label="Step size"
      lazy-rules
      :rules="[(val) => (val && val > 0) || 'Please type something']"
    >
      <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
        Simulation step size
      </q-tooltip>
    </q-input>
    <q-input
      filled
      type="number"
      v-model="stepSize"
      label="Step size"
      lazy-rules
      :rules="[(val) => (val && val > 0) || 'Please type something']"
    >
      <q-tooltip class="text-body1" :offset="[10, 10]" max-width="500px">
        Simulation step size
      </q-tooltip>
    </q-input>
    <q-checkbox
      v-model="replaceNonStandardResidues"
      label="Replace non standard residues"
    />
    <q-checkbox v-model="addMissingAtoms" label="Add missing atoms" />
    <q-checkbox v-model="addMissingHydrogens" label="Add missing hydrogens" />
    <q-checkbox v-model="frictionCoeff" label="Friction coefficient" />
    <q-checkbox v-model="ignoreMissingAtoms" label="Ignore missing atoms" />
    <q-select
      v-model="integrator"
      :options="IntegratorsValues"
      label="Integrator"
    ></q-select>
    <q-file
      filled
      bottom-slots
      accept=".pdb"
      :rules="[pdbFileRule]"
      v-model="pdbFile"
      label=".pdb file"
      counter
    >
      <template v-slot:prepend>
        <q-icon name="cloud_upload" @click.stop.prevent />
      </template>
      <template v-slot:append>
        <q-icon
          name="close"
          @click.stop.prevent="pdbFile = null"
          class="cursor-pointer"
        />
      </template>

      <template v-slot:hint> Upload .pdb file with protein structure </template>
    </q-file>
    <div>
      <q-btn
        label="Run inference"
        size="large"
        type="submit"
        color="positive"
      />
    </div>
  </q-form>
</template>
