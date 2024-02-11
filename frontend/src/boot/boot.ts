import { boot } from 'quasar/wrappers'
import {QSpinnerOrbit} from "quasar";

// "async" is optional;
// more info on params: https://v2.quasar.dev/quasar-cli/boot-files
export default boot(async ({app, router}) => {
    app.config.globalProperties.$loaderShow = (message: string) => {
        app.config.globalProperties.$q.loading.show({
            spinner: QSpinnerOrbit,
            message: message
        });
    }

    app.config.globalProperties.$loaderHide = (message: string) => {
        app.config.globalProperties.$q.loading.hide();
    }
})
