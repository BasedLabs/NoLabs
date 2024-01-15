<script lang="ts">
import { defineComponent, ref } from 'vue';
import { useSolubilityStore } from './storage';

export default defineComponent({
    setup() {
        const store = useSolubilityStore();
        const aminoAcid: string = ref(null);

        const onSubmit = async () => {
            await store.fetchInference(null, 'asdasd',aminoAcid, []);
            store.aminoAcids
        };

        const onReset = () => {
            // Your reset logic here
        };

        return {
            store,
            aminoAcid,
            onSubmit,
            onReset
        };
    }
});
</script>

<template>
    <div class="q-pa-md" style="max-width: 400px">
        <q-form @submit="onSubmit" @reset="onReset" class="q-gutter-md">
            <q-input filled v-model="aminoAcid" label="Amino acid" hint="Amino acid sequence" lazy-rules
                :rules="[val => val && val.length > 0 || 'Please type something']" />
            <div>
                <q-btn label="Submit" type="submit" color="primary" />
                <q-btn label="Reset" type="reset" color="primary" flat class="q-ml-sm" />
            </div>
        </q-form>
    </div>
</template>