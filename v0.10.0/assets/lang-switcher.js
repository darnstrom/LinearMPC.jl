(function () {
    var STORAGE_KEY = 'linearMPC-preferred-language';

    function setLanguage(switcher, lang) {
        switcher.querySelectorAll('.lang-switcher-tab').forEach(function (tab) {
            tab.classList.toggle('active', tab.dataset.lang === lang);
        });
        switcher.querySelectorAll('.lang-switcher-content').forEach(function (content) {
            content.classList.toggle('active', content.dataset.lang === lang);
        });
    }

    function switchLang(button) {
        var lang = button.dataset.lang;
        var switcher = button.closest('.lang-switcher');
        setLanguage(switcher, lang);
        try {
            localStorage.setItem(STORAGE_KEY, lang);
        } catch (e) {}
        // Sync all other switchers on the page to the same language
        document.querySelectorAll('.lang-switcher').forEach(function (sw) {
            if (sw !== switcher) { setLanguage(sw, lang); }
        });
    }

    document.addEventListener('DOMContentLoaded', function () {
        // Restore saved language preference (default: julia)
        var savedLang = 'julia';
        try {
            savedLang = localStorage.getItem(STORAGE_KEY) || 'julia';
        } catch (e) {}

        document.querySelectorAll('.lang-switcher').forEach(function (switcher) {
            setLanguage(switcher, savedLang);
        });

        // Run syntax highlighting on code blocks inside switchers
        if (typeof hljs !== 'undefined') {
            document.querySelectorAll('.lang-switcher pre code').forEach(function (el) {
                hljs.highlightElement(el);
            });
        }

        // Attach click handlers to all tab buttons
        document.querySelectorAll('.lang-switcher-tab').forEach(function (tab) {
            tab.addEventListener('click', function () {
                switchLang(this);
            });
        });
    });
})();
