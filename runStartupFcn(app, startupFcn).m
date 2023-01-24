function runStartupFcn(app, startfcn)
            ams = appdesigner.internal.service.AppManagementService.instance();
            ams.runStartupFcn(app, startfcn);
end